from __future__ import annotations

from mimosa.types import ComparisonResult

from motifhorde.comparison import UniversalMotifComparator
from motifhorde.models import GenericModel


def motif(name: str) -> GenericModel:
    return GenericModel("test", name, None, 1, {})


def test_universal_comparator_uses_mimosa_api_for_adjusted_pvalues(monkeypatch):
    calls = []

    def fake_compare_one_to_many(**kwargs):
        calls.append(kwargs)
        return [
            ComparisonResult(
                query="query",
                target="target",
                score=0.7,
                offset=0,
                orientation="++",
                metric="co",
                adj_p_value=0.03,
            )
        ]

    monkeypatch.setattr(
        "motifhorde.comparison.mimosa_compare_one_to_many",
        fake_compare_one_to_many,
    )
    comparator = UniversalMotifComparator(
        metric="co",
        comparison_criterion="p-value",
        pvalue=True,
        null_distribution="profile-null.joblib",
    )
    query = motif("query")
    target = motif("target")

    frame = comparator.compare([query], [target], sequences="sequences")

    assert calls[0]["query"] is query
    assert calls[0]["targets"] == [target]
    assert calls[0]["strategy"] == "profile"
    assert calls[0]["sequences"] == "sequences"
    assert calls[0]["comparator"] is comparator.config
    assert frame.loc[0, "adj.p-value"] == 0.03
