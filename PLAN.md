# Рефакторинг MotifHORDE 0.2.0 на Python и актуальную Mimosa

## Целевая архитектура

Ветка `refactoring` сейчас совпадает с `main`. Julia-ветка используется как
архитектурный ориентир, но код остаётся Python-проектом.

Базовая версия —
[Mimosa v1.4.6](https://github.com/ubercomrade/mimosa/tree/v1.4.6): в
`pyproject.toml` указать `mimosa-tool>=1.4.6,<1.5`, а в `uv.lock`
зафиксировать конкретную версию.

| Сейчас | После рефакторинга |
|---|---|
| `GenericModel`, изменяемый `config` | Неизменяемые `mimosa.MotifModel` |
| Собственные обёртки сканирования и FASTA | Публичный API Mimosa |
| Joblib/pickle-модели | Portable model bundles Mimosa |
| Pandas DataFrame внутри pipeline | Списки и небольшие внутренние dataclass-записи |
| Встроенный SiteGA, C++, JAR | Только внешние процессы |
| Старый Tomtom-подобный comparator | Только `mimosa.compare_many` |
| Параметры bootstrap кодируются в имени модели | Явные `task_id`, `split`, `params`, `metrics` |
| Поиск результатов SiteGA по именам файлов | Строгий версионный JSON manifest |

## Этапы реализации

### 1. Зафиксировать поведение перед миграцией

- Сохранить небольшие эталонные результаты текущих ROC/PR/pAUC,
  bootstrap-разбиения и выбора мотивов.
- Зафиксировать четыре позиционных CLI-аргумента, названия
  discovery-инструментов и детерминированную формулу seed:
  `base_seed + task_index + 1`.
- Не сохранять совместимость с Python API `GenericModel` и `.pkl`.

### 2. Сделать пакет чистым Python-проектом

- Обновить диапазоны до Python `>=3.12,<3.14`, NumPy `>=2.4.6,<3` и Mimosa
  `>=1.4.6,<1.5`.
- Удалить прямые зависимости `pandas`, `scipy`, `joblib`, `numba` и
  `pybind11`, если после миграции они нигде напрямую не импортируются.
- Удалить `setup.py`, встроенный `src/sitega`, каталог с JAR и package-data
  для бинарных файлов.
- Оставить существующий setuptools backend; не менять его без функциональной
  необходимости.
- Удалить неиспользуемые `cache.py`, `execute.py` и старый модуль реэкспорта
  `models.py`. MotifHORDE не должен дублировать API Mimosa.

### 3. Перевести доменную часть на Mimosa 1.4

- Использовать напрямую `MotifModel`, `EncodedSequences`, `read_fasta`,
  `read_model`, `write_model`, `scan`, `compare_many`, `reconstruct_pfm`,
  `TopFractionHits` и `BestPerSequence`.
- Хранить имена FASTA рядом с `EncodedSequences`; подвыборки создавать через
  `EncodedSequences.from_rows`, сохраняя порядок имён.
- Не изменять модели. Для присвоения имени встроенной модели использовать
  `dataclasses.replace`; статистику и параметры хранить отдельно.
- Ввести только две внутренние записи: `EvaluatedMotif` для модели с контекстом
  bootstrap и `MatchedPair` для результата сравнения. Удалить привязку
  статистики к строковому имени модели.
- Оставить в локальном I/O только запись FASTA/Jstacs FASTA и объединённого
  MEME PFM — этого API в Mimosa нет.

### 4. Перестроить evaluation, comparison и pipeline

- Получать позиционные оценки через `scan(..., strands="best")`.
- Перенести старую семантику ROC, PR и partial AUC в небольшие приватные
  NumPy-функции; результаты сверить с зафиксированными эталонами. При
  отсутствии валидных положительных или отрицательных оценок завершаться с
  понятной ошибкой.
- Не записывать статистику в модель.
- Для PFM использовать `reconstruct_pfm` с
  `TopFractionHits(0.10, BestPerSequence())`; не создавать собственный
  model-specific PFM-код.
- Удалить `TomtomComparator`, `GeneralMotifComparator` и
  DataFrame-преобразования. Вызывать `compare_many` для каждой query-модели с
  `total_threads=jobs`, `inner_threads=1`.
- Поддерживать только `co`, `dice`, `cosine`. Greedy-фильтрацию и дедупликацию
  выполнять стабильной сортировкой обычных списков.
- Для p-value загружать directory bundle через
  `mimosa.io.read_null_bundle`, создавать `NullDistribution` и вызывать
  `annotate_results`. BH-коррекцию применять отдельно для набора targets каждой
  query-модели.
- Сохранить отдельные рабочие каталоги bootstrap-задач, стабильный порядок
  результатов и отсутствие вложенного распараллеливания discovery-процессов.

### 5. Заменить persistence и выходной контракт

- Записывать каждую модель через `mimosa.write_model` в каталог без
  расширения: `output/<tool>/bootstrap/models/...` и
  `output/<tool>/motifs/models/...`.
- Сохранить `all_motifs_in_pfm_form.meme`.
- Сделать `statistics.json` версионным: `schema_version`, `tool`, упорядоченный
  список `motifs`; для bootstrap хранить `id`, `model`, `name`, `split`,
  `params`, `metrics`, для финальных моделей — `rank`, `model`, `name`,
  `params`, `metrics`.
- Не добавлять загрузчик pickle. Старые результаты должны быть пересозданы или
  экспортированы в MEME.
- Не перезаписывать существующие model bundle-каталоги: при конфликте
  завершаться до записи и оставлять старые данные нетронутыми.

## Внешние инструменты и CLI

- Все процессы запускать аргументным списком через `subprocess`, без shell.
  Ненулевой код должен сообщать команду, stdout и stderr.
- Приоритет разрешения executable: CLI-флаг → переменная окружения → `PATH`.
  Никаких путей внутрь пакета.
- `environment.yml` использовать как conda-first окружение приложения: сейчас
  добавить доступные `meme` и `bammmotif2`; остальные пакеты добавлять после
  публикации в conda.
- Для Dimont и Slim временно оставить
  `--dimont-jar`/`HORDEMOTIFS_DIMONT_JAR` и
  `--slim-jar`/`HORDEMOTIFS_SLIM_JAR`. Встроенного fallback не будет.
- Добавить `--sitega-command`/`HORDEMOTIFS_SITEGA_COMMAND`; SiteGA всегда
  запускается как процесс.
- Conda-рецепты SiteGA, Dimont и Slim должны жить в их проектах или
  feedstock-репозиториях, не в MotifHORDE.

CLI сохраняет четыре позиционных аргумента и инструменты `streme`, `meme`,
`bamm`, `dimont`, `slim`, `sitega`, но:

- удалить `-c`/`--comparator`, `--tomtom-metric` и `--pfm-mode`;
- сохранить `--mimosa-metric`, но ограничить значениями `co`, `dice`,
  `cosine`;
- `--mimosa-null-distribution` теперь принимает каталог bundle, а не
  joblib-файл;
- исправить SiteGA-параметры: `--lpd` означает максимальную LPD-дистанцию
  `1..6`, а новый `--features` задаёт сетку количества признаков. Старое
  использование `--lpd` как числа признаков не поддерживать.

Внешний SiteGA должен принять `--manifest <path>` и записать:

```json
{
  "schema_version": 1,
  "status": "success",
  "message": "",
  "models": [
    {
      "model_file": "foreground_mat1.mat",
      "model_type": "sitega",
      "length": 12,
      "name": "Sitega-1"
    }
  ]
}
```

MotifHORDE проверяет schema/status, относительность пути и отсутствие выхода из
output-каталога, существование файла, тип и длину через Mimosa, сохраняет
manifest-порядок и возвращает не больше запрошенного числа моделей. Поиск
`.mat` по шаблону удалить.

## Проверка и критерии готовности

- Unit-тесты Mimosa API: FASTA, подвыборки, scan, все поддерживаемые model
  readers, bundle round-trip и PFM reconstruction.
- Регрессионные тесты ROC/PR/pAUC и выбора финальных мотивов.
- Тесты comparison для трёх метрик, score/p-value порогов, стабильного порядка
  и совместимости null bundle.
- Fake-process тесты каждого discovery-адаптера, включая аргументы,
  env/CLI precedence, ненулевой exit code и отсутствующий output.
- Отдельные manifest-тесты SiteGA: неверная версия, failure status, path
  traversal, отсутствующая модель, неверные type/length и сохранение порядка.
- Повторный запуск bootstrap с одним seed должен давать одинаковые task IDs,
  порядок и JSON.
- Реальные discovery-тесты оставить opt-in; обычная CI не должна требовать
  внешних программ.
- Финальные команды: `uv lock`, `uv sync --locked`, `uv run pytest`,
  `uv run ruff check`, `uv run ruff format --check`, `uv build`,
  `git diff --check`.
- Проверить содержимое wheel: никаких `.jar`, `.cpp`, `.so`, pickle/joblib и
  SiteGA implementation.

## Допущения

- Это breaking release `0.2.0`; compatibility facade не создаётся.
- Julia-код не переносится и не включается в Python-пакет.
- ChipMunk в этом рефакторинге не добавляется: его conda-пакет можно
  использовать при отдельном добавлении нового discovery-адаптера.
- До появления всех conda-пакетов явные command/JAR-флаги остаются
  поддерживаемым переходным механизмом.
