"""Synthetic decision-gate tests; no real tag recommendation or functional data."""
from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("tag_vector_tools", ROOT / "scripts/vector_tools.py")
vt = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(vt)


def setup(tmp_path):
    (tmp_path / "source.fa").write_text(">synthetic_test\nATGAAATAA\n")
    data = {
        "schema_version": "1.0", "construct_id": "tag_gate_test",
        "topology": "linear", "replication_class": "de_novo",
        "sources": [{
            "id": "source", "file": "source.fa", "format": "fasta", "evidence": {
                "source_url": None, "accession_version": None,
                "boundary_basis": "synthetic software test, not biological material",
                "evidence_class": "test_fixture", "reviewed": True,
            },
        }],
        "assembly": [{"source_id": "source", "start": 1, "end": 9, "strand": 1}],
        "features": [{
            "id": "target", "type": "CDS", "segments": [[1, 9]], "strand": 1,
            "expected_protein": "MK", "complete_cds": True,
        }],
        "new_expression_target_ids": ["target"],
        "unresolved_issues": [],
    }
    return tmp_path / "construct.json", data


def no_tag(target="target"):
    return {
        "target_feature_id": target, "design_context": "new_design", "choice": "no_tag",
        "decision_source": "user_explicit", "user_response": "Unit-test synthetic explicit no-tag choice",
        "purpose": None, "selected_tag": None, "fusion_site": "none",
        "literature_evidence_ids": [], "sequence_source": None, "tag_feature_ids": [],
        "plan_confirmed": False,
    }


def add_tag():
    decision = no_tag()
    decision.update({
        "choice": "add_tag", "user_response": "Unit-test synthetic explicit add-tag choice",
        "purpose": "test protocol only", "selected_tag": "SYNTHETIC_NOT_A_REAL_TAG",
        "fusion_site": "C_terminal", "literature_evidence_ids": ["TEST_EVIDENCE_NOT_A_PAPER"],
        "sequence_source": "test_fixture/source.fa", "tag_feature_ids": ["tag_block"],
        "plan_confirmed": True,
        "fusion_plan": c_terminal_fusion_plan(),
    })
    return decision


def build(tmp_path, data):
    path = tmp_path / "construct.json"
    vt.dump_json(path, data)
    return vt.build(path, tmp_path / "out", True)


def c_terminal_fusion_plan():
    return {
        "assembly_mode": "continuous_orf",
        "components": [
            {"feature_id": "gene_body", "role": "target", "source_id": "gene_native",
             "source_segments": [[1, 9]], "source_strand": 1,
             "stop_action": "remove_for_fusion", "expected_peptide": "MK"},
            {"feature_id": "flexible_linker", "role": "flexible_linker", "source_id": "linker_source",
             "source_segments": [[1, 15]], "source_strand": 1,
             "stop_action": "none_present", "expected_peptide": "GGGGS"},
            {"feature_id": "tag_block", "role": "tag", "source_id": "tag_source",
             "source_segments": [[1, 6]], "source_strand": 1,
             "stop_action": "retained_as_final", "expected_peptide": "A"},
        ],
        "linker_flexibility_basis": "Synthetic GS-rich protocol fixture; no functional claim.",
        "linker_evidence_ids": ["TEST_EVIDENCE_NOT_A_PAPER"],
        "final_stop": {
            "codon": "TAA", "provider": "last_component", "source_id": "tag_source",
            "source_segments": [[4, 6]], "source_strand": 1,
        },
    }


def add_fixture_source(data, tmp_path, source_id, sequence):
    (tmp_path / (source_id + ".fa")).write_text(f">{source_id}\n{sequence}\n")
    source = copy.deepcopy(data["sources"][0])
    source["id"] = source_id
    source["file"] = source_id + ".fa"
    data["sources"].append(source)


def insert_test_tag_feature(data, tmp_path):
    sequence = "ATGAAA" + "GGTGGTGGTGGTTCT" + "GCT" + "TAA"
    (tmp_path / "source.fa").write_text(">synthetic_test\n" + sequence + "\n")
    data["assembly"][0]["end"] = len(sequence)
    data["features"][0]["segments"] = [[1, len(sequence)]]
    data["features"][0]["expected_protein"] = "MKGGGGSA"
    for fid, start, end in [("gene_body", 1, 6), ("flexible_linker", 7, 21), ("tag_block", 22, 24)]:
        data["features"].append({
            "id": fid, "type": "misc_feature", "segments": [[start, end]], "strand": 1,
            "qualifiers": {"label": "SYNTHETIC FUSION FIXTURE " + fid},
        })
    add_fixture_source(data, tmp_path, "gene_native", "ATGAAATAA")
    add_fixture_source(data, tmp_path, "linker_source", "GGTGGTGGTGGTTCT")
    add_fixture_source(data, tmp_path, "tag_source", "GCTTAA")


def test_de_novo_without_decision_is_blocked(tmp_path):
    _, data = setup(tmp_path)
    with pytest.raises(Exception, match="tag_decisions"):
        build(tmp_path, data)
    assert not (tmp_path / "out/construct.gb").exists()


def test_no_tag_requires_explicit_user_choice_but_not_tag_search(tmp_path):
    _, data = setup(tmp_path)
    data["tag_decisions"] = [no_tag()]
    qc = build(tmp_path, data)
    assert qc["tag_decision_checks"][0]["choice"] == "no_tag"


def test_pending_does_not_default_to_no_tag(tmp_path):
    _, data = setup(tmp_path)
    decision = no_tag()
    decision.update(choice="pending", decision_source="not_asked", user_response=None)
    data["tag_decisions"] = [decision]
    with pytest.raises(ValueError, match="pending"):
        build(tmp_path, data)


def test_reported_tag_cannot_replace_user_decision_for_de_novo(tmp_path):
    _, data = setup(tmp_path)
    decision = no_tag()
    decision.update(design_context="replicated_original", decision_source="reported_design", user_response=None)
    data["tag_decisions"] = [decision]
    with pytest.raises(ValueError, match="explicit user"):
        build(tmp_path, data)


@pytest.mark.parametrize("field,value,error", [
    ("plan_confirmed", False, "not confirmed"),
    ("literature_evidence_ids", [], "literature evidence"),
    ("sequence_source", None, "sequence_source"),
    ("selected_tag", None, "selected_tag"),
    ("purpose", None, "purpose"),
    ("tag_feature_ids", [], "tag feature IDs"),
    ("fusion_site", "pending", "fusion site"),
])
def test_added_tag_requires_confirmed_evidence_backed_plan(tmp_path, field, value, error):
    _, data = setup(tmp_path)
    insert_test_tag_feature(data, tmp_path)
    decision = add_tag()
    decision[field] = value
    data["tag_decisions"] = [decision]
    with pytest.raises(ValueError, match=error):
        build(tmp_path, data)


def test_confirmed_synthetic_plan_preserves_complete_translation(tmp_path):
    _, data = setup(tmp_path)
    insert_test_tag_feature(data, tmp_path)
    data["tag_decisions"] = [add_tag()]
    qc = build(tmp_path, data)
    assert qc["tag_decision_checks"][0]["status"] == "declared_plan_consistency_pass"
    assert qc["translation_checks"][0]["translation"] == "MKGGGGSA"
    assert qc["status"] == "test_fixture_only"


def test_each_de_novo_cds_requires_its_own_decision(tmp_path):
    _, data = setup(tmp_path)
    extra = copy.deepcopy(data["features"][0])
    extra["id"] = "second_target"
    data["features"].append(extra)
    data["new_expression_target_ids"].append("second_target")
    data["tag_decisions"] = [no_tag()]
    with pytest.raises(ValueError, match="second_target"):
        build(tmp_path, data)


def test_partial_redesign_new_product_also_needs_real_user_decision(tmp_path):
    _, data = setup(tmp_path)
    data["replication_class"] = "partial_redesign"
    decision = no_tag()
    decision["decision_source"] = "not_asked"
    data["tag_decisions"] = [decision]
    with pytest.raises(ValueError, match="explicit user"):
        build(tmp_path, data)


def test_existing_replication_without_new_tag_plan_remains_compatible(tmp_path):
    _, data = setup(tmp_path)
    data["replication_class"] = "literature_reconstruction"
    del data["new_expression_target_ids"]
    qc = build(tmp_path, data)
    assert qc["tag_decision_checks"] == []


def test_tag_annotation_must_be_within_target_cds_same_strand(tmp_path):
    _, data = setup(tmp_path)
    insert_test_tag_feature(data, tmp_path)
    data["features"][-1]["strand"] = -1
    data["tag_decisions"] = [add_tag()]
    with pytest.raises(ValueError, match="same strand"):
        build(tmp_path, data)


def test_no_tag_metadata_cannot_silently_select_tag(tmp_path):
    _, data = setup(tmp_path)
    decision = no_tag()
    decision["selected_tag"] = "UNCONFIRMED"
    data["tag_decisions"] = [decision]
    with pytest.raises(ValueError, match="no_tag conflicts"):
        build(tmp_path, data)


def test_duplicate_or_unknown_target_is_blocked(tmp_path):
    _, data = setup(tmp_path)
    data["tag_decisions"] = [no_tag(), no_tag()]
    with pytest.raises(ValueError, match="duplicate"):
        build(tmp_path, data)
    data["tag_decisions"] = [no_tag("unknown")]
    with pytest.raises(ValueError, match="actual CDS"):
        build(tmp_path, data)


def test_from_scratch_gate_does_not_infer_an_unannotated_target(tmp_path):
    _, data = setup(tmp_path)
    data["features"] = []
    data["tag_decisions"] = [no_tag()]
    with pytest.raises(ValueError, match="actual CDS"):
        build(tmp_path, data)


def test_tag_workflow_is_conditional_and_evidence_specific():
    hub = (ROOT / "SKILL.md").read_text()
    rules = (ROOT / "references/tag-selection.md").read_text()
    report = (ROOT / "templates/design-report.md").read_text()
    assert "references/tag-selection.md" in hub
    assert "用户需要标签后才启动" in rules
    assert "具体支持实验" in rules
    assert "不能预设" in rules and "通用最优" in rules
    assert "标签决策与文献推荐" in report
    assert "不添加" in report and "不将模型偏好" in report


def test_example_is_pending_and_does_not_make_user_choice():
    data = json.loads((ROOT / "templates/tag-decision.example.json").read_text())
    assert data["choice"] == "pending"
    assert data["decision_source"] == "not_asked"
    assert data["selected_tag"] is None and not data["plan_confirmed"]


def test_backbone_cds_is_not_a_default_tag_question_target(tmp_path):
    _, data = setup(tmp_path)
    marker = copy.deepcopy(data["features"][0])
    marker["id"] = "synthetic_backbone_marker_fixture"
    data["features"].append(marker)
    data["tag_decisions"] = [no_tag()]
    qc = build(tmp_path, data)
    assert len(qc["tag_decision_checks"]) == 1


def test_tag_component_cds_does_not_recursively_require_a_tag_decision(tmp_path):
    _, data = setup(tmp_path)
    insert_test_tag_feature(data, tmp_path)
    tag = next(f for f in data["features"] if f["id"] == "tag_block")
    tag.update({"type": "CDS", "expected_protein": "A", "complete_cds": False})
    data["tag_decisions"] = [add_tag()]
    qc = build(tmp_path, data)
    assert len(qc["tag_decision_checks"]) == 1
