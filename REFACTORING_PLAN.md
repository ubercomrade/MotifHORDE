# План исправлений по результатам ревью кода

Дата: 2026-09-21. Ветка: `refactoring`. База: коммит `bd5905a` (rm plans).

Исходное ревью: Ruff и тесты чистые (42 passed / 3 skipped), CLI работает.
Проблемы сгруппированы в 8 этапов ниже. Каждый этап — отдельный коммит; после
каждого этапа тесты должны быть зелёными без изменения ожиданий существующих
тестов (кроме явно указанных адаптаций).

---

## Принятые решения по спорным моментам

Спорные вопросы были заданы, ответов не поступило — приняты рекомендуемые
варианты. Любой из них можно переиграть до начала соответствующего этапа.

| # | Вопрос | Принято | Альтернатива |
|---|--------|---------|--------------|
| 1 | Судьба ~15 МБ неиспользуемых тестовых данных | Удалить всё неиспользуемое; история recoverable из git | Удалить, но сохранить `PEAKS*.fa` для будущих интеграционных тестов |
| 2 | Фиктивные `nsites=20` в MEME-экспорте | Вычислять реальное значение из PFM | Оставить константу, добавить комментарий |
| 3 | Объём рефакторинга `pipeline.py` | Умеренный: `ComparisonSettings` + дедупликация отбора, структуру модулей не менять | Минимальный (только гигиена) или глубокий (разбить на selection/persistence) |
| 4 | Валидация null-бандла при каждом `compare_records` | Убрать из `compare_records`; валидировать один раз в пайплайне (с sequences — в `run()`) | Оставить как есть (корректность важнее скорости) |

---

## Этап 0. Гигиена git-репозитория

Решение №1. Изменений кода нет — только индекс git.

Удалить из индекса (данные остаются в истории):

- `tests/__pycache__/` — 13 файлов `.pyc`, включая байткод удалённых тестов
  старого API (`test_tomtom_comparator`, `test_correlation_comparator`,
  `test_motif_comparator`). Правила `__pycache__/` (строка 73) и `*.pyc`
  (строка 156) уже есть в `.gitignore` — достаточно `git rm -r --cached`.
- Неиспользуемые тестовые данные (ни один тест и документ на них не
  ссылается; единственный потребитель `test_data` —
  `tests/test_external_fullrun.py:36`, использует только `small_pipeline/`):
  - `tests/test_data/PEAKS035260_GATA2_O09100_MACS2.fa` и `_gb.fa`
  - `tests/test_data/PEAKS035588_AR_P19091_MACS2.fa` и `_gb.fa`
  - `tests/test_data/gimme/` (логи и intermediate стороннего инструмента)
  - `tests/test_data/results/` (сгенерированные результаты streme,
    включая `motifs.meme`)
  - `tests/test_data/hist_pass`, `tests/test_data/prc_pass`
  - `tests/test_data/run.sh`, `tests/test_data/gimme.sh`
  - `train.fa` (корень репозитория)
  - `tests/test_output.fa`

Убрать личные пути из `tests/test_external_fullrun.py:54-69`: дефолты
`/home/anon/Programs/Jstacs/Dimont.jar` и `.../SlimDimont.jar` заменить на
обязательные переменные окружения (поведение `_require_jar`: без env — skip).
Попутно убрать захардкоженный путь из `id=` параметризации.

Приёмка:

```bash
git ls-files | grep -c '\.pyc$'          # → 0
git ls-files | grep -v small_pipeline | grep test_data   # → пусто
uv run pytest                            # 42 passed, 3 skipped
```

---

## Этап 1. Мёртвый код и остатки рефакторинга

Файл `src/motifhorde/discovery.py`:

1. `_build_dimont_args` (строки 369–399) не вызывается: `DimontDiscoveryTool.
   discover` строит аргументы через `_build_jstacs_args(..., extra=[ess=...])`
   напрямую. Восстановить симметрию со Slim: `discover` вызывает
   `_build_dimont_args`, прямой вызов `_build_jstacs_args` с `extra=` убираем.
   Дублирования знания «ess-флаг Dimont» больше нет.
2. `StremeDiscoveryTool.__init__(nmotifs=5, ...)` (строки 180–183): параметр
   `nmotifs` и атрибут `self.nmotifs` нигде не читаются — `discover()` берёт
   `number_of_motifs` из аргументов. Удалить параметр и атрибут.
   Сопутствующее изменение: `cli.py:240`
   `StremeDiscoveryTool(args.nmotifs, args.streme_command)` →
   `StremeDiscoveryTool(args.streme_command)`.
   Тесты конструктор не ломаются (используют `command=` ключевым аргументом).

Файл `src/motifhorde/pipeline.py`:

3. `import tempfile` внутри `_temporary_discovery_dir` (строка 415) → в шапку
   файла.
4. `from .comparison import comparison_value` внутри `_candidate_value`
   (строка 504) → убрать; `comparison_value` добавить в существующий импорт
   из `.comparison` в шапке.

Приёмка: `grep -rn "_build_dimont_args\|self\.nmotifs" src/` — только вызов и
определение; `grep -n "    import \|    from " src/motifhorde/*.py` — пусто;
ruff + pytest зелёные.

---

## Этап 2. ComparisonSettings и дедупликация жадного отбора

Решение №3 (умеренный объём). Публичный конструктор `DeNovoPipeline` и CLI
не меняются.

### 2.1. Frozen-dataclass настроек сравнения

Новый `ComparisonSettings` в `comparison.py`:

```python
@dataclass(frozen=True, slots=True)
class ComparisonSettings:
    metric: str                  # mimosa-метрика: co | dice | cosine
    criterion: str               # score | p-value
    threshold: float
    null_distribution: NullDistribution | None = None
    search_range: int = 10
    jobs: int = 1
```

- Валидация в `__post_init__`: метрика из `SUPPORTED_METRICS`, критерий
  `score`/`p-value`, `jobs >= 1`, `search_range >= 0`, порог конечен;
  `p-value` требует null (перенос проверок из `DeNovoPipeline.__init__`
  строк 62–97 и `compare_records` строк 113–118 — валидация перестаёт
  дублироваться в трёх местах).
- Метод `compare(self, query, targets, sequences) -> list[MatchedPair]` —
  тонкая обёртка над `compare_records`, чтобы вызывающий код не протаскивает
  5 аргументов.
- `DeNovoPipeline` сохраняет текущие параметры конструктора, но внутри
  хранит `self.comparison: ComparisonSettings`.

### 2.2. Единая запись финального мотива

Вместо нетипизированных `dict[str, Any]` в `_select_final_motifs` /
`_deduplicate_final_motifs` / `_save_results` — frozen-dataclass:

```python
@dataclass(frozen=True, slots=True)
class FinalMotif:
    model: MotifModel
    name: str
    params: MappingProxyType   # нормализация как в EvaluatedMotif
    metrics: MappingProxyType
```

### 2.3. Один жадный отбор вместо двух

`_select_nonredundant_motifs` (pipeline.py:471, работает с `EvaluatedMotif`)
и `_deduplicate_final_motifs` (pipeline.py:511, работает с dict) — идентичный
алгоритм «сортируй по метрике ↓, добавляй если не похож ни на один уже
отобранный». Заменить одной функцией:

```python
def _select_nonredundant(
    items: Sequence[tuple[MotifModel, float, T]],   # (модель, ключ сортировки, payload)
    sequences: EncodedSequences,
    settings: ComparisonSettings,
) -> list[T]
```

Оба call-site упаковывают свои записи в кортежи и распаковывают результат.
`_compare_bootstrap_motifs` заодно теряет два одинаковых 9-аргументных блока
(строки 239–260) — заменяются на цикл `for split_group in (odd, even)` или
два вызова с settings.

### 2.4. Мелочи, попутно

- `_average_metric` (pipeline.py:467) не должен пересобирать dict всех четырёх
  метрик на каждое сравнение в `sorted()` (строка 314): вычислять нужную
  метрику напрямую.
- `_candidate_value` упраздняется: используется `comparison_value` из шапки
  (после этапа 1) — при необходимости оставить локальную lambda.

Приёмка: тесты `test_pipeline_output` (порядок результатов, независимость
от jobs, seed) проходят **без изменения ожиданий**. Добавить юнит-тест
единой функции отбора: два мотива выше порога сходства → остаётся лучший.

Риск этапа: средний (логика отбора). Защита — стабильные ожидания
существующих тестов порядка.

---

## Этап 3. Таблица инструментов в cli.py

`check_dependencies` (cli.py:160–234) — 75 строк if/elif, дублирующих знание
«инструмент → команда → env-переменная», которое уже есть в адаптерах
`discovery.py`. Заменить декларативной таблицей:

```python
@dataclass(frozen=True)
class CommandRequirement:
    label: str
    default: str          # DEFAULT_*_COMMAND из external.py
    cli_attr: str         # имя атрибуса Namespace
    env_var: str

TOOL_COMMANDS: dict[str, tuple[CommandRequirement, ...]] = {
    "streme": (CommandRequirement("STREME", DEFAULT_STREME_COMMAND, "streme_command", "HORDEMOTIFS_STREME_COMMAND"),),
    "meme":   (...),
    "bamm":   (STREME + BaMMmotif),
    "sitega": (...),
}
# dimont/slim: java-команда + jar через resolve_existing_path (отдельные записи)
```

- `check_dependencies` и `setup_discovery_tool` обходят таблицу; if/elif-цепь
  и `_dependency_error` (тривиальный враппер над `raise RuntimeError`,
  cli.py:156) удаляются.
- Тексты ошибок сохранить дословно (на них могут завязываться пользователи).
- `_resolve_jobs` (cli.py:132): расставить скобки
  `return (os.cpu_count() or 1) if jobs == -1 else jobs` — поведение то же.

Приёмка: новые тесты `check_dependencies` (monkeypatch env/PATH) на все шесть
инструментов и ветку отсутствующей зависимости; `test_cli.py` продолжает
проходить.

---

## Этап 4. Валидация null-бандла один раз

Решение №4. Сейчас `compare_records` (comparison.py:120–126) валидирует бандл
при каждом вызове, включая пересчёт `normalization_fingerprint` и
`sequence_fingerprint` во вложенных циклах отбора и `_select_final_motifs`.

- Убрать вызов `validate_null_distribution` из `compare_records`.
- `DeNovoPipeline.__init__` валидирует бандл без sequences (как сейчас,
  строки 92–97); дополнительно в `run()` — одна валидация с `sequences=peaks`
  перед первым сравнением (сохраняет привязку бандла к последовательностям).
- `validate_null_distribution` остаётся публичной для прямых пользователей.

Адаптация тестов: `test_comparison.py:119–145` (несовместимый search_range,
metric и т.д. через `compare_records`) — перенести ожидания `pytest.raises`
на явный вызов `validate_null_distribution`. Тест конструктора пайплайна с
несовместимым бандлом (`test_pipeline_output.py:128+`) остаётся как есть.

Документация: заметка в `MIGRATION.md` — `compare_records` больше не
перевалидирует бандл на каждый вызов; прямой вызов с чужим бандлом должен
валидироваться явно.

Риск: прямой вызов `compare_records` с несовместимым бандлом перестаёт
защищаться. Приемлемо: функция внутренняя, пайплайн валидирует на входе.

---

## Этап 5. Честные метаданные MEME + мелкие дефекты

Решение №2.

### 5.1. Реальный nsites в `write_meme` (io.py:80)

PFM реконструируется через `TopFractionHits(0.10, BestPerSequence())` —
каждая последовательность даёт не более одного хита, все колонки PFM
суммируются в одно и то же число сайтов. Изменения:

- `write_meme(pfms, records, path)`, где `record = (name, length, nsites)`;
  в заголовке `nsites= {int}` вместо константы 20. `E= 0` оставить
  (E-value неизвестен, это допустимый placeholder формата MEME).
- `pipeline._save_results` вычисляет
  `nsites = int(round(float(pfm.sum(axis=0).max())))` для каждого мотива.
- Тест: записанный файл содержит согласованный nsites (равен сумме колонок);
  защита от расползания колонок (если суммы различаются — ошибка валидации
  в `write_meme`).

### 5.2. Прочее

- `external.py:56–59`: сообщение об ошибке запуска процесса содержит пустой
  блок `stdout:`. Заменить на
  `f"Command failed to start: {shlex.join(args)}\nerror: {exc}"` — stdout/
  stderr у незапустившегося процесса не существует.
- `discovery.py:586–591`: `raise ValueError` внутри `try` перехватывается
  собственным `except ValueError` и перебрасывается. Оставить один raise;
  `try` сузить до вызова `os.path.commonpath` (его `ValueError` на Windows —
  отдельный случай, сейчас замаскирован).

Приёмка: тест `write_meme` на nsites; ruff; pytest.

---

## Этап 6. Типизация и аннотации

Без изменения поведения, только точность контрактов:

- `discovery.py`: аннотировать `-> list[MotifModel]` у `discover()` в ABC и
  всех шести адаптерах; в docstring ABC перечислить ожидаемые kwargs
  (`length` — везде, `order` — bamm, `lpd`/`features`/`seed` — sitega).
  Сигнатуру `*args, **kwargs` не трогаем (решение №3 — умеренный объём).
- `pipeline.py`: `null_distribution: NullDistribution | str | os.PathLike |
  None = None` в конструкторе; аннотации у `_ensure_pfm -> np.ndarray`,
  `_temporary_discovery_dir -> tempfile.TemporaryDirectory`,
  `_iter_param_grid -> Iterator[dict[str, Any]]`.
- `comparison.py`: аннотировать `sequences: EncodedSequences` в
  `compare_records`; поля `MatchedPair` оставить `Any` (принимают и
  `MotifModel`, и `EvaluatedMotif`, и `FinalMotif`), описать в docstring.
- `cli.py`: `setup_discovery_tool(args) -> MotifDiscoveryTool`.

Приёмка: ruff зелёный, `uv build` собирается, поведение не изменилось
(pytest без правок).

---

## Этап 7. Покрытие тестами CLI-слоя

Сейчас `tests/test_cli.py` — один тест парсера. Добавить:

- `parse_range`: одиночное значение; список через запятую; диапазон со шагом
  (`"8-20-4"` → `[8, 12, 16, 20]`); `step=0` → ValueError; пустой результат
  → ValueError.
- `validate_args`: `--jobs 0` / `--jobs -2` → SystemExit; `p-value` без
  `--mimosa-null-distribution` → SystemExit; `--lpd 7` для sitega →
  SystemExit; дефолт порога заполняется (после этапа 2 — из
  `ComparisonSettings`, убрать дублирующее заполнение в `validate_args`
  cli.py:142, оставить в одном месте).
- `check_dependencies` (после этапа 3): каждый инструмент с доступной
  командой проходит; недоступная команда → RuntimeError с именем метки.
- `main_cli --help` → SystemExit(0).

Приёмка: `uv run pytest` — новых скипов нет, всё детерминированно и быстро
(fake-процессы, реальных вызовов нет).

---

## Этап 8. Финальная верификация и документация

1. `README.md`:
   - раздел про MEME-экспорт: nsites теперь реальный;
   - таблица env-переменных без изменений (проверить актуальность после
     этапа 3).
2. `MIGRATION.md`: заметка про `compare_records` и null-валидацию (этап 4).
3. Полный прогон по чек-листу AGENTS.md:

```bash
uv run ruff format --check .
uv run ruff check .
uv run pytest
uv run motifhorde --help
uv build
git diff --check
git ls-files | grep -Ec '\.pyc$|train\.fa|test_output'   # → 0
```

4. Проверить содержимое wheel (`uv build` + unzip -l): нет тестовых данных,
   нет egg-info-мусора, только `src/motifhorde`.

---

## Порядок выполнения и зависимости

```text
Этап 0 (гигиена)          — независим, делать первым (меньше шума в диффах)
Этап 1 (мёртвый код)      — независим
Этап 2 (ComparisonSettings) — этап 3 частично опирается на него
Этап 3 (таблица CLI)      — после 1 (Streme-конструктор)
Этап 4 (null-валидация)   — после 2 (settings); меняет тесты comparison
Этап 5 (MEME + мелочи)    — независим
Этап 6 (типизация)        — после 2 (аннотации финальных структур)
Этап 7 (тесты CLI)        — после 3 (тестируем таблицу)
Этап 8 (верификация)      — последним
```

Оценка рисков: этапы 0, 1, 5, 6 — низкий; 3, 7 — низкий/средний; 2, 4 —
средний (логика отбора и контракт сравнения), закрываются стабильными
ожиданиями существующих тестов и явными адаптациями, перечисленными выше.
