"""Sequence-observation and design-review contract tests, not functional prediction."""
import csv
import importlib.util
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("sequence_review", ROOT / "scripts/sequence_review.py")
sr = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(sr)


def feature(kind, start, end, strand=1, label="test"):
    return SeqFeature(SimpleLocation(start, end, strand=strand), type=kind, qualifiers={"label": [label]})


def test_atg_candidate_is_observed_not_functionally_called():
    rec = SeqRecord(Seq("CCCATGAAATAACCC"), id="leader_fixture")
    rec.features = [feature("5'UTR", 0, len(rec))]
    result = sr.review(rec)
    candidate = result["leader_review"][0]["candidates"][0]
    assert candidate["leader_start_nt"] == 4
    assert candidate["length_aa_before_stop"] == 2
    assert candidate["translation_or_function_validated"] is False
    assert result["source_modified"] is False


def test_unresolved_stop_does_not_extend_past_annotation():
    rec = SeqRecord(Seq("CCCATGAAACCC"), id="partial_leader")
    rec.features = [feature("5'UTR", 0, 9)]
    item = sr.review(rec)["leader_review"][0]["candidates"][0]
    assert item["stop_leader_start_nt"] is None
    assert "not_resolved" in item["candidate_class"]


def test_negative_strand_respects_declared_orientation():
    coding = "CCCATGAAATAACCC"
    rec = SeqRecord(Seq(coding).reverse_complement(), id="negative_leader")
    rec.features = [feature("5'UTR", 0, len(rec), -1)]
    item = sr.review(rec)["leader_review"][0]["candidates"][0]
    assert item["start_triplet_positions_1based"] == [12, 11, 10]


def test_exact_repeats_do_not_imply_a_biological_defect():
    part = "ACGTACGTACGTAACCGGTT"
    rec = SeqRecord(Seq(part + "AAAA" + part), id="repeat_fixture")
    rec.features = [
        feature("promoter", 0, len(part), label="p1"),
        feature("promoter", len(part) + 4, len(rec), label="p2"),
    ]
    result = sr.review(rec)
    assert len(result["exact_regulatory_repeats"]) == 1
    assert "no silencing/recombination outcome inferred" in result["flags"][0]["observation"]


def test_double_annotation_of_same_locus_is_not_a_repeat():
    rec = SeqRecord(Seq("ACGT" * 8), id="annotation_aliases")
    rec.features = [feature("promoter", 0, 24, label="a"), feature("promoter", 0, 24, label="b")]
    assert sr.review(rec)["exact_regulatory_repeats"] == []


def test_leader_review_is_not_lost_under_an_overlapping_promoter_annotation():
    rec = SeqRecord(Seq("CCCATGAAATAACCC"), id="overlapping_annotation")
    rec.features = [feature("promoter", 0, len(rec)), feature("5'UTR", 0, len(rec))]
    assert len(sr.review(rec)["leader_review"][0]["candidates"]) == 1


def test_missing_leader_is_not_guessed_from_promoter():
    rec = SeqRecord(Seq("ATGAAATAA"), id="missing_leader")
    rec.features = [feature("promoter", 0, len(rec))]
    result = sr.review(rec)
    assert not result["leader_review"]
    assert any(f["category"] == "not_assessed" and f["feature"] == "5prime_leader" for f in result["flags"])


def test_non_aug_is_explicitly_out_of_scope():
    rec = SeqRecord(Seq("CCCACGAAATAACCC"), id="non_AUG_fixture")
    rec.features = [feature("5'UTR", 0, len(rec))]
    result = sr.review(rec)
    assert result["leader_review"][0]["candidates"] == []
    assert "non-AUG" in result["limitations"][0]


def test_ambiguous_sequence_is_not_reported_as_clean():
    with pytest.raises(ValueError, match="ambiguous"):
        sr.review(SeqRecord(Seq("ACGTNN"), id="ambiguous"))


def test_fuzzy_and_unknown_directions_are_not_silently_resolved():
    rec = SeqRecord(Seq("CCCATGAAATAACCC"), id="fuzzy")
    rec.features = [SeqFeature(SimpleLocation(BeforePosition(0), len(rec), strand=1), type="5'UTR")]
    with pytest.raises(ValueError, match="fuzzy"):
        sr.review(rec)
    rec.features = [feature("5'UTR", 0, len(rec), None)]
    assert any(x["category"] == "unassessed_annotation" for x in sr.review(rec)["flags"])


def test_generate_preserves_source_and_outputs_facts_only(tmp_path):
    rec = SeqRecord(Seq("CCCATGAAATAACCC"), id="source_fixture")
    rec.annotations["molecule_type"] = "DNA"
    rec.features = [feature("5'UTR", 0, len(rec))]
    source = tmp_path / "source.gb"
    SeqIO.write(rec, source, "genbank")
    original = source.read_bytes()
    result = sr.generate(source, tmp_path / "review")
    assert source.read_bytes() == original
    assert result["status"] == "sequence_observations_only"
    assert (tmp_path / "review/review_flags.tsv").exists()


def test_design_review_has_all_twelve_topics_and_no_default_pass():
    with (ROOT / "templates/design-review.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [r["review_id"] for r in rows] == [f"R{i:02d}" for i in range(1, 13)]
    assert all(r["review_status"] == "gap" for r in rows)


def test_annotated_review_does_not_auto_redesign_parts():
    rules = (ROOT / "references/evidence-driven-review.md").read_text()
    assert "不自动删除" in rules and "不修改任何序列" in rules
    assert "unresolved_issues" in rules and "2A" in rules
    assert "单株外推" in rules and "启动子–leader–CDS–终止区域" in rules
