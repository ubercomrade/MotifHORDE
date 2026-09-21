# План закрытия остаточных проблем рефакторинга

Этот документ дополняет `PLAN.md`. `PLAN.md` остаётся архитектурным контрактом
релиза 0.2.0, а здесь перечислены подтверждённые аудитом незакрытые работы и
порядок их приёмки.

## 1. Исправить p-value comparison

Файлы: `src/motifhorde/comparison.py`, `tests/test_comparison.py`.

- Создавать `NullDistribution(**read_null_bundle(path))`, а не передавать сырой
  словарь в `annotate_results`.
- До аннотации проверять совместимость bundle с метрикой и параметрами profile
  comparison.
- Использовать только публичный API Mimosa. Если Mimosa 1.4 не предоставляет
  достаточной публичной проверки совместимости, сначала добавить её в Mimosa,
  выпустить совместимую версию 1.4.x и затем поднять минимальную версию
  `mimosa-tool`; не копировать fingerprinting и не импортировать приватный API.
- Добавить тесты успешной загрузки настоящего directory bundle, score- и
  p-value-порогов, BH-коррекции отдельно для targets каждой query,
  несовместимых параметров и повреждённого bundle.

Критерий готовности: полный pipeline проходит с
`--comparison-criterion p-value`.

## 2. Удалить локальный MEME-парсер

Файлы: `src/motifhorde/discovery.py`, тесты discovery-адаптеров.

- Удалить `_meme_stdout_for_mimosa` и связанные parsing helpers.
- Получать от MEME стандартный файл, который напрямую принимает
  `mimosa.read_model`.
- Если Mimosa не читает стандартный MEME output, исправить reader в Mimosa,
  выпустить совместимую версию 1.4.x и обновить нижнюю границу зависимости.
  Не переносить parser обратно в MotifHORDE.
- Добавить fake-process тест с несколькими мотивами и проверкой порядка, имён и
  длины.

Критерий готовности: MotifHORDE не разбирает модельные матрицы; весь model I/O
проходит через Mimosa.

## 3. Сделать SiteGA manifest строгим

Файлы: `src/motifhorde/discovery.py`, `tests/test_sitega_manifest.py`.

- Сначала валидировать все записи `models`, затем применять лимит
  `number_of_motifs`.
- Проверять наличие и типы всех обязательных полей, включая `message`.
- Сохранить проверки относительного POSIX-пути, containment и symlink.
- Добавить тесты неверной версии, failure status, пустого/невалидного списка,
  отсутствующей модели, неверных `model_type` и length, абсолютного пути,
  traversal, symlink, ошибки после лимита, manifest order и ограничения числа
  возвращаемых моделей.

## 4. Восстановить регрессионные тесты поведения

Файлы: `tests/test_evaluation_bootstrap.py`,
`tests/test_pipeline_selection.py`.

- Зафиксировать точные эталонные `auROC`, `auPRC`, `pauROC` и `pauPRC` для
  небольших score-векторов, включая ties и крайние FPR.
- Проверить ошибки при отсутствии валидных положительных или отрицательных
  оценок.
- Восстановить тесты стабильной сортировки, inclusive thresholds, локальности
  сопоставления по параметрам, выбора full-data мотива по обеим bootstrap-моделям
  и глобальной дедупликации.
- Выполнить два bootstrap-запуска в разных каталогах с одним seed и сравнить
  task IDs, порядок records и `statistics.json` байт-в-байт.
- Проверить одинаковый результат при `jobs=1` и `jobs=2`.

Эталонные значения должны происходить из дорефакторингового поведения, а не
вычисляться тестируемыми функциями во время теста.

## 5. Восстановить fake-process coverage адаптеров

Предпочтительно использовать один параметризованный
`tests/test_discovery_adapters.py`, без отдельного файла на каждый инструмент.

Для `streme`, `meme`, `bamm`, `dimont`, `slim` и `sitega` проверить:

- точный argument list и запуск без shell;
- приоритет CLI -> environment -> `PATH`;
- отсутствие вложенного распараллеливания: один внутренний thread;
- ненулевой exit code с командой, кодом, stdout и stderr в ошибке;
- отсутствующий или пустой output;
- неверную модель и несовместимую длину;
- стабильный порядок результатов.

Отдельно проверить Jstacs FASTA и разрешение JAR через CLI/environment.

## 6. Закрыть Mimosa и integration contracts

- Добавить компактные fixtures для readers MEME/PWM, BaMM, Dimont, Slim и
  SiteGA. Для каждого непосредственно вызвать `mimosa.read_model` и проверить
  тип и длину.
- Сохранить проверки FASTA, подвыборок, scan, bundle round-trip и PFM
  reconstruction.
- Восстановить opt-in real-tool smoke tests с маркерами `external` и `fullrun`.
- Обновить ожидания integration-тестов: directory bundles вместо `.pkl`,
  внешний SiteGA manifest и отсутствие bundled JAR.
- Обычный `uv run pytest` не должен запускать внешние программы.

## 7. Документация и финальная приёмка

- Добавить `MIGRATION.md` с удалёнными Python API, отказом от `.pkl`, новыми
  bundle paths, null directory bundle и контрактом внешнего SiteGA.
- Обновить README и отметить выполненные пункты в этом документе только после
  прохождения соответствующих тестов.
- Выполнить:

```bash
uv lock
uv sync --locked
uv run pytest
uv run ruff check .
uv run ruff format --check .
uv run motifhorde --help
uv build
git diff --check
```

- Проверить wheel и sdist на отсутствие `.jar`, `.cpp`, `.so`, pickle/joblib,
  SiteGA implementation и сгенерированных model bundles.

## Порядок выполнения

1. p-value path.
2. Граница MEME/Mimosa.
3. SiteGA manifest.
4. Регрессии evaluation и selection.
5. Fake-process и opt-in integration tests.
6. Документация и финальная приёмка.

Новые зависимости, compatibility wrappers и дополнительные архитектурные слои
для этих работ не требуются.

## Статус приёмки

- [x] p-value comparison: directory bundle, public `NullDistribution`,
  compatibility checks, thresholds, and per-query BH correction.
- [x] MEME boundary: standard output is passed directly to Mimosa; local
  matrix parsing is removed.
- [x] strict SiteGA manifest validation and limit-after-validation behavior.
- [x] ROC/PR/pAUC references, invalid-class errors, stable bootstrap ordering,
  seeded reproducibility, and `jobs=1`/`jobs=2` parity.
- [x] fake-process coverage for discovery adapters, command resolution,
  Jstacs FASTA/JAR handling, and Mimosa reader contracts.
- [x] migration documentation and packaging checks.
- [ ] real external-tool smoke tests: implemented as `external`/`fullrun`
  tests, not run without `HORDEMOTIFS_RUN_FULLRUN=1` and installed tools.
