# Стратегия рефакторинга: Переход на Numba-совместимый RaggedScores

## 1. Архитектурные изменения
Текущая реализация `RaggedScores` в `hordemotifs/functions.py` использует массив `values` формы `(n_seq, Lmax)`, что приводит к значительным затратам памяти из-за padding.

Новое решение в `hordemotifs/ragged.py` основано на:
- **Упакованном представлении (CSR-like):** Единый одномерный массив `data` + массив смещений `offsets`.
- **Numba jitclass:** Класс `RaggedScoresNumba` доступен напрямую внутри `@njit` функций.

## 2. Пошаговый план миграции

### Шаг 1: `hordemotifs/models.py`
- Добавить импорт `RaggedScoresNumba` и `ragged_from_list` из `.ragged`.
- В методе `get_scores` добавить возможность возврата нового формата.
- В методе `_calculate_threshold_table` использовать `ragged.data` напрямую (больше нет нужды в масках).

### Шаг 2: `hordemotifs/functions.py`
- Заменить старую функцию `ragged` на новую реализацию без `np.full((n, Lmax), ...)`.
- Переписать `_fast_overlap_kernel_masked` и `_fast_cj_kernel_masked` для работы с `RaggedScoresNumba`. Новые версии ядер избавятся от 4-х лишних параметров (`len1_max`, `len2_max`, `valid1`, `valid2`).

### Шаг 3: `hordemotifs/io.py`
- Оптимизировать `_read_sequences` так, чтобы последовательности сразу записывались в плоский буфер, полностью исключая стадию создания `List[np.ndarray]`.

### Шаг 4: `hordemotifs/pipeline.py`
- Обновить пайплайны для передачи объектов `RaggedScoresNumba` в функции обработки.

## 3. Таблица изменений

| Файл | Точка интеграции | Тип изменения |
| --- | --- | --- |
| `models.py` | `get_scores` | Возврат `RaggedScoresNumba` |
| `functions.py` | `_fast_overlap_kernel` | Удаление масок и 2D индексов |
| `ragged.py` | — | Новый модуль с `jitclass` |
| `io.py` | `read_fasta` | Прямое заполнение плоского массива |

## 4. Пример "До и После"

### До (2D массивы + Маски)
```python
for i in range(n_seq):
    row_ptr1 = i * len1_max
    for j in range(overlap):
        v1 = s1_flat[row_ptr1 + j] # Сложный индекс
```

### После (jitclass)
```python
for i in range(ragged.num_sequences):
    s = ragged.get_slice(i)
    for j in range(len(s)):
        v = s[j] # Прямой доступ к срезу