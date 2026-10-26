"""Synthetic ORF/junction tests only; no actual plant tag or linker optimization."""
import copy
import json

import pytest
from Bio.Seq import Seq

from test_tag_selection import setup, add_tag, insert_test_tag_feature, add_fixture_source, build


LINKER_DNA = "GGTGGTGGTGGTTCT"


def prepared_c(tmp_path):
    _, data = setup(tmp_path)
    insert_test_tag_feature(data, tmp_path)
    data["tag_decisions"] = [add_tag()]
    return data


def feature(data, name):
    return next(f for f in data["features"] if f["id"] == name)


def test_c_terminal_removes_gene_stop_and_retains_final_tag_stop(tmp_path):
    data = prepared_c(tmp_path)
    qc = build(tmp_path, data)
    fusion = qc["tag_decision_checks"][0]["fusion_junctions"]
    assert fusion["status"] == "continuous_orf_and_junctions_pass"
    assert fusion["components"][0]["source_terminal_stop"] == "TAA"
    assert fusion["components"][0]["stop_action"] == "remove_for_fusion"
    assert fusion["components"][1]["encoded_peptide"] == "GGGGS"
    assert fusion["components"][1]["linker_dna"] == LINKER_DNA
    assert len(fusion["junctions"]) == 2
    assert fusion["components"][-1]["stop_action"] == "retained_as_final"
    assert "flexible_linker" in (tmp_path / "out/expression_map.svg").read_text()


def test_missing_fusion_plan_blocks_new_tag_finalization(tmp_path):
    data = prepared_c(tmp_path)
    del data["tag_decisions"][0]["fusion_plan"]
    with pytest.raises(ValueError, match="fusion_plan"):
        build(tmp_path, data)
    assert not (tmp_path / "out/construct.gb").exists()


def test_upstream_stop_cannot_be_retained(tmp_path):
    data = prepared_c(tmp_path)
    data["tag_decisions"][0]["fusion_plan"]["components"][0]["stop_action"] = "retained_as_final"
    with pytest.raises(ValueError, match="upstream terminal stop"):
        build(tmp_path, data)


def test_stop_removal_must_not_be_claimed_for_stopless_source(tmp_path):
    data = prepared_c(tmp_path)
    (tmp_path / "gene_native.fa").write_text(">gene_native\nATGAAA\n")
    data["tag_decisions"][0]["fusion_plan"]["components"][0]["source_segments"] = [[1, 6]]
    with pytest.raises(ValueError, match="no source terminal stop"):
        build(tmp_path, data)


def test_stopless_source_is_recorded_as_none_present(tmp_path):
    data = prepared_c(tmp_path)
    (tmp_path / "gene_native.fa").write_text(">gene_native\nATGAAA\n")
    part = data["tag_decisions"][0]["fusion_plan"]["components"][0]
    part.update(source_segments=[[1, 6]], stop_action="none_present")
    qc = build(tmp_path, data)
    assert qc["tag_decision_checks"][0]["fusion_junctions"]["components"][0]["source_terminal_stop"] is None


def test_linker_must_be_declared_flexible_at_target_tag_boundary(tmp_path):
    data = prepared_c(tmp_path)
    data["tag_decisions"][0]["fusion_plan"]["components"][1]["role"] = "other"
    with pytest.raises(ValueError, match="flexible linker is required"):
        build(tmp_path, data)


def test_linker_peptide_must_match_actual_dna(tmp_path):
    data = prepared_c(tmp_path)
    data["tag_decisions"][0]["fusion_plan"]["components"][1]["expected_peptide"] = "GGGGA"
    with pytest.raises(ValueError, match="peptide does not match"):
        build(tmp_path, data)


def test_linker_codon_boundary_is_checked(tmp_path):
    data = prepared_c(tmp_path)
    feature(data, "flexible_linker")["segments"] = [[7, 20]]
    with pytest.raises(ValueError, match="codon boundary"):
        build(tmp_path, data)


def test_internal_stop_in_linker_is_rejected(tmp_path):
    data = prepared_c(tmp_path)
    bad = "GGTTAAGGTGGTTCT"
    (tmp_path / "source.fa").write_text(">synthetic_test\nATGAAA" + bad + "GCTTAA\n")
    (tmp_path / "linker_source.fa").write_text(">linker_source\n" + bad + "\n")
    with pytest.raises(ValueError, match="internal stop"):
        build(tmp_path, data)


def test_final_stop_requires_actual_terminal_codon(tmp_path):
    data = prepared_c(tmp_path)
    (tmp_path / "source.fa").write_text(">synthetic_test\nATGAAA" + LINKER_DNA + "GCTCAA\n")
    with pytest.raises(ValueError, match="final stop codon"):
        build(tmp_path, data)


def test_declared_stop_must_match_source_and_orf(tmp_path):
    data = prepared_c(tmp_path)
    data["tag_decisions"][0]["fusion_plan"]["final_stop"]["codon"] = "TAG"
    with pytest.raises(ValueError, match="does not match"):
        build(tmp_path, data)


def test_declared_n_terminal_must_match_component_order(tmp_path):
    data = prepared_c(tmp_path)
    data["tag_decisions"][0]["fusion_site"] = "N_terminal"
    with pytest.raises(ValueError, match="N-to-C component order"):
        build(tmp_path, data)


def test_components_cannot_hide_unannotated_scar_or_gap(tmp_path):
    data = prepared_c(tmp_path)
    data["tag_decisions"][0]["fusion_plan"]["components"] = list(reversed(
        data["tag_decisions"][0]["fusion_plan"]["components"]
    ))
    with pytest.raises(ValueError, match="order, coverage"):
        build(tmp_path, data)


def test_vector_can_supply_final_stop_if_tag_source_has_none(tmp_path):
    data = prepared_c(tmp_path)
    (tmp_path / "tag_source.fa").write_text(">tag_source\nGCT\n")
    plan = data["tag_decisions"][0]["fusion_plan"]
    plan["components"][-1].update(source_segments=[[1, 3]], stop_action="none_present")
    add_fixture_source(data, tmp_path, "vector_stop", "TAA")
    plan["final_stop"].update(provider="vector", source_id="vector_stop", source_segments=[[1, 3]])
    qc = build(tmp_path, data)
    assert qc["tag_decision_checks"][0]["fusion_junctions"]["final_stop"]["provider"] == "vector"


def test_last_component_provider_requires_retained_source_stop(tmp_path):
    data = prepared_c(tmp_path)
    data["tag_decisions"][0]["fusion_plan"]["components"][-1]["stop_action"] = "remove_for_fusion"
    with pytest.raises(ValueError, match="not recorded as retained_as_final"):
        build(tmp_path, data)


def test_negative_strand_fusion_uses_coding_orientation(tmp_path):
    data = prepared_c(tmp_path)
    full = "ATGAAA" + LINKER_DNA + "GCTTAA"
    (tmp_path / "source.fa").write_text(">synthetic_test\n" + str(Seq(full).reverse_complement()) + "\n")
    n = len(full)
    for f in data["features"]:
        f["segments"] = [[n - end + 1, n - start + 1] for start, end in f["segments"]]
        f["strand"] = -1
    assert build(tmp_path, data)["translation_checks"][0]["translation"] == "MKGGGGSA"


def test_negative_source_cds_stop_is_read_on_declared_strand(tmp_path):
    data = prepared_c(tmp_path)
    (tmp_path / "gene_native.fa").write_text(">gene_native\nTTATTTCAT\n")
    data["tag_decisions"][0]["fusion_plan"]["components"][0]["source_strand"] = -1
    assert build(tmp_path, data)["tag_decision_checks"][0]["fusion_junctions"]["components"][0]["source_terminal_stop"] == "TAA"


def prepared_n(tmp_path):
    data = prepared_c(tmp_path)
    full = "ATGGCT" + LINKER_DNA + "ATGAAATAA"
    (tmp_path / "source.fa").write_text(">synthetic_test\n" + full + "\n")
    (tmp_path / "tag_source.fa").write_text(">tag_source\nATGGCTTAA\n")
    data["assembly"][0]["end"] = len(full)
    data["features"][0].update(segments=[[1, len(full)]], expected_protein="MAGGGGSMK")
    feature(data, "tag_block")["segments"] = [[1, 6]]
    feature(data, "flexible_linker")["segments"] = [[7, 21]]
    feature(data, "gene_body")["segments"] = [[22, 27]]
    decision = data["tag_decisions"][0]
    decision["fusion_site"] = "N_terminal"
    plan = decision["fusion_plan"]
    gene, linker, tag = plan["components"]
    tag.update(source_segments=[[1, 9]], expected_peptide="MA", stop_action="remove_for_fusion")
    gene["stop_action"] = "retained_as_final"
    plan["components"] = [tag, linker, gene]
    plan["final_stop"].update(source_id="gene_native", source_segments=[[7, 9]])
    return data


def test_n_terminal_removes_tag_stop_but_can_retain_target_stop(tmp_path):
    data = prepared_n(tmp_path)
    qc = build(tmp_path, data)
    fusion = qc["tag_decision_checks"][0]["fusion_junctions"]
    assert fusion["components"][0]["role"] == "tag"
    assert fusion["components"][0]["stop_action"] == "remove_for_fusion"
    assert fusion["components"][-1]["role"] == "target"
    assert fusion["components"][-1]["stop_action"] == "retained_as_final"
    assert qc["translation_checks"][0]["translation"] == "MAGGGGSMK"


def test_n_terminal_does_not_keep_upstream_tag_stop(tmp_path):
    data = prepared_n(tmp_path)
    data["tag_decisions"][0]["fusion_plan"]["components"][0]["stop_action"] = "retained_as_final"
    with pytest.raises(ValueError, match="upstream terminal stop"):
        build(tmp_path, data)


def prepared_internal(tmp_path):
    data = prepared_c(tmp_path)
    full = "ATGAAA" + LINKER_DNA + "GTT" + LINKER_DNA + "GCTTAA"
    (tmp_path / "source.fa").write_text(">synthetic_test\n" + full + "\n")
    (tmp_path / "gene_native.fa").write_text(">gene_native\nATGAAAGCTTAA\n")
    (tmp_path / "tag_source.fa").write_text(">tag_source\nGTTTAA\n")
    data["assembly"][0]["end"] = len(full)
    data["features"] = [{
        "id": "target", "type": "CDS", "segments": [[1, len(full)]], "strand": 1,
        "expected_protein": "MKGGGGSVGGGGSA", "complete_cds": True,
    }] + [
        {"id": fid, "type": "misc_feature", "segments": [[start, end]], "strand": 1}
        for fid, start, end in [
            ("gene_n", 1, 6), ("linker_n", 7, 21), ("tag_block", 22, 24),
            ("linker_c", 25, 39), ("gene_c", 40, 42),
        ]
    ]
    decision = data["tag_decisions"][0]
    decision["fusion_site"] = "internal"
    base = decision["fusion_plan"]["components"]
    gene_n = copy.deepcopy(base[0])
    gene_n.update(feature_id="gene_n", source_segments=[[1, 6]], stop_action="none_present")
    linker_n = copy.deepcopy(base[1])
    linker_n["feature_id"] = "linker_n"
    tag = copy.deepcopy(base[2])
    tag.update(stop_action="remove_for_fusion", expected_peptide="V")
    linker_c = copy.deepcopy(linker_n)
    linker_c["feature_id"] = "linker_c"
    gene_c = copy.deepcopy(gene_n)
    gene_c.update(feature_id="gene_c", source_segments=[[7, 12]], stop_action="retained_as_final", expected_peptide="A")
    decision["fusion_plan"]["components"] = [gene_n, linker_n, tag, linker_c, gene_c]
    decision["fusion_plan"]["final_stop"].update(source_id="gene_native", source_segments=[[10, 12]])
    return data


def test_internal_tag_checks_both_flexible_junctions(tmp_path):
    data = prepared_internal(tmp_path)
    qc = build(tmp_path, data)
    assert qc["translation_checks"][0]["translation"] == "MKGGGGSVGGGGSA"


def test_internal_tag_cannot_omit_one_flexible_junction(tmp_path):
    data = prepared_internal(tmp_path)
    data["tag_decisions"][0]["fusion_plan"]["components"][3]["role"] = "other"
    with pytest.raises(ValueError, match="every target/tag junction"):
        build(tmp_path, data)


def test_existing_replicated_tag_is_not_silently_reengineered(tmp_path):
    data = prepared_c(tmp_path)
    data["replication_class"] = "literature_reconstruction"
    del data["new_expression_target_ids"]
    decision = data["tag_decisions"][0]
    decision.update(design_context="replicated_original", decision_source="reported_design", user_response=None)
    del decision["fusion_plan"]
    qc = build(tmp_path, data)
    assert qc["tag_decision_checks"][0]["fusion_junctions"] == "not_audited_legacy_replication"


def test_fusion_template_is_explicitly_incomplete():
    from test_tag_selection import ROOT
    template = json.loads((ROOT / "templates/fusion-plan.example.json").read_text())
    assert template["components"] == []
    assert template["final_stop"]["codon"] is None


def test_allowed_alternative_start_is_interpreted_as_initiator_only_at_orf_start(tmp_path):
    data = prepared_c(tmp_path)
    (tmp_path / "source.fa").write_text(">synthetic_test\nTTGAAA" + LINKER_DNA + "GCTTAA\n")
    (tmp_path / "gene_native.fa").write_text(">gene_native\nTTGAAATAA\n")
    qc = build(tmp_path, data)
    assert qc["translation_checks"][0]["translation"] == "MKGGGGSA"


@pytest.mark.parametrize("stop", ["TAA", "TAG", "TGA"])
def test_upstream_source_stop_is_identified_from_the_actual_sequence(tmp_path, stop):
    data = prepared_c(tmp_path)
    (tmp_path / "gene_native.fa").write_text(">gene_native\nATGAAA" + stop + "\n")
    qc = build(tmp_path, data)
    assert qc["tag_decision_checks"][0]["fusion_junctions"]["components"][0]["source_terminal_stop"] == stop


def test_fusion_can_be_assembled_directly_from_stop_trimmed_source_parts(tmp_path):
    data = prepared_c(tmp_path)
    data["assembly"] = [
        {"source_id": "gene_native", "start": 1, "end": 6, "strand": 1},
        {"source_id": "linker_source", "start": 1, "end": 15, "strand": 1},
        {"source_id": "tag_source", "start": 1, "end": 6, "strand": 1},
    ]
    assert build(tmp_path, data)["translation_checks"][0]["translation"] == "MKGGGGSA"
