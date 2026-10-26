"""Synthetic software tests only. No biological validation or actual Evo2 inference."""
from __future__ import annotations

import copy
import csv
import importlib.util
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


vt = load("vector_tools")
ej = load("evo2_job")


def cds(segments=None, strand=1):
    return {"id": "test_cds", "type": "CDS", "segments": segments or [[1, 9]],
            "strand": strand, "expected_protein": "MK", "complete_cds": True}


def setup_manifest(tmp_path, sequence="ATGAAATAA", features=None):
    (tmp_path / "source.fa").write_text(">test_source\n" + sequence + "\n")
    data = {
        "schema_version": "1.0", "construct_id": "software_test", "topology": "linear",
        "replication_class": "literature_reconstruction",
        "sources": [{"id": "source", "file": "source.fa", "format": "fasta", "evidence": {
            "source_url": None, "accession_version": None, "boundary_basis": "Synthetic software test",
            "evidence_class": "test_fixture", "reviewed": True}}],
        "assembly": [{"source_id": "source", "start": 1, "end": len(sequence), "strand": 1}],
        "features": [cds()] if features is None else features,
        "protected_regions": [], "restriction_sites": [], "required_feature_ids": [],
        "unresolved_issues": [],
    }
    path = tmp_path / "manifest.json"
    vt.dump_json(path, data)
    return path, data


def test_full_export_roundtrip_and_svg(tmp_path):
    path, _ = setup_manifest(tmp_path)
    qc = vt.build(path, tmp_path / "out", True)
    assert qc["status"] == "test_fixture_only"
    assert qc["file_roundtrip"] == "pass"
    assert qc["translation_checks"][0]["translation"] == "MK"
    record = vt.read_record(tmp_path / "out/construct.gb", "genbank")
    assert str(record.seq) == "ATGAAATAA"
    assert len(record.features) == 1
    for name in ("plasmid_map.svg", "expression_map.svg"):
        assert ET.parse(tmp_path / "out" / name).getroot().tag.endswith("svg")
    assert (tmp_path / "out/features.tsv").exists()


@pytest.mark.parametrize("sequence", ["ATGNAATAA", "AUGAAAUAA", ""])
def test_unambiguous_dna_required(tmp_path, sequence):
    path, _ = setup_manifest(tmp_path, sequence, features=[])
    with pytest.raises(Exception):
        vt.build(path, tmp_path / "out", True)


def test_fixture_not_real_design(tmp_path):
    path, _ = setup_manifest(tmp_path)
    with pytest.raises(ValueError, match="Test fixture"):
        vt.build(path, tmp_path / "out")
    assert not (tmp_path / "out/construct.gb").exists()


def test_missing_source_and_unresolved_issue_block(tmp_path):
    path, data = setup_manifest(tmp_path)
    data["unresolved_issues"] = ["missing promoter boundary"]
    vt.dump_json(path, data)
    with pytest.raises(ValueError, match="Unresolved"):
        vt.build(path, tmp_path / "out", True)


def test_hash_mismatch(tmp_path):
    path, data = setup_manifest(tmp_path)
    data["sources"][0]["sequence_sha256"] = "0" * 64
    vt.dump_json(path, data)
    with pytest.raises(ValueError, match="hash mismatch"):
        vt.build(path, tmp_path / "out", True)


@pytest.mark.parametrize(
    "sequence,feature,topology",
    [
        ("ATGAAAGTCCAGTAA", cds([[1, 6], [13, 15]]), "linear"),
        ("TTATTTCAT", cds([[1, 9]], -1), "linear"),
        ("TAAGGGGATGAAA", cds([[8, 13], [1, 3]]), "circular"),
    ],
)
def test_join_negative_and_cross_origin_cds(tmp_path, sequence, feature, topology):
    path, data = setup_manifest(tmp_path, sequence, [feature])
    data["topology"] = topology
    vt.dump_json(path, data)
    qc = vt.build(path, tmp_path / "out", True)
    assert qc["translation_checks"][0]["translation"] == "MK"


def test_wrong_protein_blocks_export(tmp_path):
    path, data = setup_manifest(tmp_path)
    data["features"][0]["expected_protein"] = "ML"
    vt.dump_json(path, data)
    with pytest.raises(ValueError, match="protein mismatch"):
        vt.build(path, tmp_path / "out", True)


def test_coordinates_and_duplicate_features(tmp_path):
    path, data = setup_manifest(tmp_path)
    data["features"][0]["segments"] = [[1, 10]]
    vt.dump_json(path, data)
    with pytest.raises(ValueError, match="coordinates"):
        vt.build(path, tmp_path / "out", True)
    data["features"] = [cds(), cds()]
    vt.dump_json(path, data)
    with pytest.raises(ValueError, match="Duplicate"):
        vt.build(path, tmp_path / "out", True)


def test_exact_plasmid_reference_verification(tmp_path):
    path, data = setup_manifest(tmp_path)
    data["replication_class"] = "exact_plasmid"
    vt.dump_json(path, data)
    with pytest.raises(ValueError, match="replication reference"):
        vt.build(path, tmp_path / "bad", True)
    data["replication_reference_source_id"] = "source"
    vt.dump_json(path, data)
    assert "pass" in vt.build(path, tmp_path / "good", True)["reference_check"]


def test_circular_equivalence_and_restriction_origin():
    assert vt.equivalent("GTAC", "ACGT", True)
    assert vt.equivalent(str(Seq("ATGAAA").reverse_complement()), "ATGAAA", True)
    assert not vt.equivalent("GTAC", "ACGT", False)
    assert any(x["crosses_origin"] for x in vt.site_hits("AAGC", "CAA", True))
    assert not vt.site_hits("AAGC", "CAA", False)


def test_protected_region_and_synonymous_edits(tmp_path):
    _, data = setup_manifest(tmp_path)
    data["protected_regions"] = [{"id": "start", "segments": [[1, 3]], "strand": 1}]
    candidate, report = vt.apply_edits(data, "ATGAAATAA", [{"pos": 6, "ref": "A", "alt": "G"}])
    assert candidate == "ATGAAGTAA" and report["status"] == "pass"
    with pytest.raises(ValueError, match="Protected region"):
        vt.apply_edits(data, "ATGAAATAA", [{"pos": 1, "ref": "A", "alt": "G"}])
    with pytest.raises(ValueError, match="Protein"):
        vt.apply_edits(data, "ATGAAATAA", [{"pos": 5, "ref": "A", "alt": "C"}])


def test_equal_length_and_ref_matching(tmp_path):
    _, data = setup_manifest(tmp_path)
    with pytest.raises(ValueError, match="REF mismatch"):
        vt.apply_edits(data, "ATGAAATAA", [{"pos": 6, "ref": "C", "alt": "G"}])
    with pytest.raises(ValueError, match="equal-length"):
        vt.apply_edits(data, "ATGAAATAA", [{"pos": 6, "ref": "A", "alt": "GG"}])
    report = vt.compare_sequences(data, "ATGAAATAA", "ATGAAAATAA")
    assert report["status"] == "fail" and "incomparable" in report["comparability"]


def test_inherited_cds_is_preserved_and_constrained(tmp_path):
    path, data = setup_manifest(tmp_path, features=[])
    record = SeqRecord(Seq("ATGAAATAA"), id="test_source")
    record.annotations["molecule_type"] = "DNA"
    record.features = [SeqFeature(SimpleLocation(0, 9, strand=1), type="CDS", qualifiers={"translation": ["MK"]})]
    SeqIO.write(record, tmp_path / "source.gb", "genbank")
    data["sources"][0].update({"file": "source.gb", "format": "genbank"})
    vt.dump_json(path, data)
    record, _, _ = vt.assemble(vt.manifest(path), tmp_path, True)
    assert record.features[0].qualifiers["translation"] == ["MK"]
    report = vt.compare_sequences(data, record.seq, "ATGACATAA", record.features)
    assert report["status"] == "fail"


def test_reverse_source_inheritance(tmp_path):
    path, data = setup_manifest(tmp_path, "TTATTTCAT")
    data["assembly"][0]["strand"] = -1
    vt.dump_json(path, data)
    assert vt.build(path, tmp_path / "out", True)["translation_checks"][0]["translation"] == "MK"


def test_partial_source_feature_must_be_reannotated(tmp_path):
    path, data = setup_manifest(tmp_path, features=[])
    record = SeqRecord(Seq("ATGAAATAA"), id="test_source")
    record.annotations["molecule_type"] = "DNA"
    record.features = [SeqFeature(SimpleLocation(0, 9, strand=1), type="CDS", qualifiers={"translation": ["MK"]})]
    SeqIO.write(record, tmp_path / "source.gb", "genbank")
    data["sources"][0].update({"file": "source.gb", "format": "genbank"})
    data["assembly"][0]["end"] = 6
    vt.dump_json(path, data)
    with pytest.raises(ValueError, match="crosses cut"):
        vt.build(path, tmp_path / "out", True)


def make_profile(tmp_path, verified=True):
    profile = {
        "profile_id": "MOCK_UNIT_TEST_NOT_EVO2", "verified": verified,
        "model_type": "evo2_7b", "model_revision": "MOCK", "code_revision": "MOCK",
        "argv_template": [sys.executable, str(ROOT / "tests/mock_backend.py"), "{input}", "{output}", "{weights_path}"],
        "env_template": {}, "device": "cuda:0", "score_type": "mean_log_likelihood",
        "score_columns": {"id": "sequence_name", "value": "probability"},
        "output_delimiter": ",", "context_limit_nt": 100,
        "limits": {"max_items": 30, "max_item_nt": 100, "timeout_seconds": 30},
        "simulation": True,
    }
    path = tmp_path / "profile.json"
    ej.dump(path, profile)
    return path, profile


def make_job(tmp_path, verified=True, candidates=None):
    profile, _ = make_profile(tmp_path, verified)
    input_file = tmp_path / "input.fa"
    input_file.write_text(candidates or ">reported_baseline\nATGAAATAAC\n>candidate\nATGAAGTAAC\n")
    weights = tmp_path / "path with spaces" / "unit-test-only.pt"
    weights.parent.mkdir()
    weights.write_text("Not a model; command protocol unit test only.")
    out = tmp_path / "job"
    job = ej.prepare(input_file, profile, str(weights), out, 6, 4)
    return out / "job.json", job, weights


def test_runtime_weight_path_not_required_during_skill_creation(tmp_path):
    profile, _ = make_profile(tmp_path)
    input_file = tmp_path / "input.fa"
    input_file.write_text(">x\nACGT\n")
    with pytest.raises(ValueError, match="Ask the user"):
        ej.prepare(input_file, profile, "", tmp_path / "job")
    # Preparation accepts a remote path without claiming it exists in this cloud.
    job = ej.prepare(input_file, profile, "/remote/hpc/later supplied.pt", tmp_path / "prepared", 4, 2)
    assert job["execution_status"] == "not_executed"


def test_explicit_windows_cover_tail(tmp_path):
    path, job, _ = make_job(tmp_path, candidates=">reported_baseline\nACGTACGTACG\n")
    covered = set()
    for item in job["items"]:
        covered.update(range(item["start"], item["end"] + 1))
    assert covered == set(range(1, 12))
    assert max(x["length"] for x in job["items"]) == 6


def test_unverified_profile_is_blocked(tmp_path):
    path, _, _ = make_job(tmp_path, False)
    with pytest.raises(ValueError, match="unverified"):
        ej.run(path)


def test_dry_run_never_executes_and_path_is_single_argv(tmp_path):
    path, _, weights = make_job(tmp_path)
    plan = ej.run(path)
    assert plan["argv"][-1] == str(weights)
    assert plan["status"] == "dry_run"
    assert not (path.parent / "run_receipt.json").exists()
    assert not (path.parent / "raw_scores.csv").exists()


def test_hash_tampering_blocks_execution(tmp_path):
    path, _, _ = make_job(tmp_path)
    (path.parent / "evo2_input.fasta").write_text(">changed\nACGT\n")
    with pytest.raises(ValueError, match="hash changed"):
        ej.run(path)


def test_command_template_needs_real_weight_binding(tmp_path):
    path, job, _ = make_job(tmp_path)
    job["profile"]["argv_template"] = ["easyevo2", "score", "{input}", "--output", "{output}"]
    ej.dump(path, job)
    with pytest.raises(ValueError, match="weights_path"):
        ej.run(path)


def test_shell_template_rejected(tmp_path):
    path, job, _ = make_job(tmp_path)
    job["profile"]["argv_template"] = ["bash", "-c", "command {input} {output} {weights_path}"]
    ej.dump(path, job)
    with pytest.raises(ValueError, match="Shell interpreter"):
        ej.run(path)


def test_simulated_backend_roundtrip_is_not_model_validation(tmp_path):
    path, _, _ = make_job(tmp_path)
    receipt = ej.run(path, True)
    assert receipt["status"] == "simulation_completed"
    with pytest.raises(ValueError, match="Simulation is not Evo2"):
        ej.collect(path, "reported_baseline")
    summary = ej.collect(path, "reported_baseline", allow_simulation=True)
    assert summary["execution_status"] == "simulation_completed"
    assert summary["calibration_status"] == "not_calibrated"
    with (path.parent / "evo2_scores.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(float(row["delta_score"]) == 0 for row in rows)


def test_receipt_required_for_collection(tmp_path):
    path, _, _ = make_job(tmp_path)
    with pytest.raises(ValueError, match="No actual run receipt"):
        ej.collect(path, "reported_baseline", True)


def test_different_parent_length_not_ranked(tmp_path):
    path, _, _ = make_job(tmp_path, candidates=">reported_baseline\nACGTACGTAC\n>candidate\nACGTACGTACG\n")
    ej.run(path, True)
    ej.collect(path, "reported_baseline", True)
    with (path.parent / "evo2_scores.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(row["delta_score"] == "NA" for row in rows if row["parent_id"] == "candidate")


def test_output_tampering_and_nonfinite_scores_rejected(tmp_path):
    path, _, _ = make_job(tmp_path)
    ej.run(path, True)
    raw = path.parent / "raw_scores.csv"
    raw.write_text(raw.read_text().replace("-1.0", "nan"))
    with pytest.raises(ValueError, match="output changed"):
        ej.collect(path, "reported_baseline", True)
    receipt_path = path.parent / "run_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["output_sha256"] = ej.sha(raw)
    ej.dump(receipt_path, receipt)
    with pytest.raises(ValueError, match="nonfinite"):
        ej.collect(path, "reported_baseline", True)


def test_resource_limit_rejects_unbounded_task(tmp_path):
    profile_path, profile = make_profile(tmp_path)
    profile["limits"]["max_items"] = 1
    ej.dump(profile_path, profile)
    file = tmp_path / "input.fa"
    file.write_text(">x\nACGTACGTACGT\n")
    with pytest.raises(ValueError, match="Item count"):
        ej.prepare(file, profile_path, "/remote/model.pt", tmp_path / "out", 4, 2)


def test_packaged_example_is_not_executable_design():
    with pytest.raises(Exception):
        vt.manifest(ROOT / "templates/construct.example.json")
    profile = json.loads((ROOT / "templates/hpc-profile.example.json").read_text())
    assert profile["weights_path"] is None and profile["argv_template"] is None and not profile["verified"]


def test_cli_failure_exports_qc_not_fake_sequence(tmp_path):
    path, data = setup_manifest(tmp_path)
    data["unresolved_issues"] = ["unknown source"]
    vt.dump_json(path, data)
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/vector_tools.py"), "build",
         "--manifest", str(path), "--out", str(tmp_path / "out"), "--allow-test-fixtures"],
        capture_output=True, text=True,
    )
    assert result.returncode == 1
    assert (tmp_path / "out/qc_report.json").exists()
    assert not (tmp_path / "out/construct.fasta").exists()


def test_cds_without_direction_is_not_validated(tmp_path):
    record = SeqRecord(Seq("ATGAAATAA"), id="ambiguous_direction")
    record.features = [SeqFeature(SimpleLocation(0, 9), type="CDS", qualifiers={"translation": ["MK"]})]
    with pytest.raises(ValueError, match="unresolved CDS direction"):
        vt.check_features(record)


def test_identical_sequence_score_inconsistency_is_rejected(tmp_path):
    path, _, _ = make_job(
        tmp_path, candidates=">reported_baseline\nACGTACGTAC\n>no_change_control\nACGTACGTAC\n"
    )
    ej.run(path, True)
    raw = path.parent / "raw_scores.csv"
    with raw.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    rows[-1]["probability"] = "-2.0"
    with raw.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sequence_name", "probability"])
        writer.writeheader()
        writer.writerows(rows)
    receipt_path = path.parent / "run_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["output_sha256"] = ej.sha(raw)
    ej.dump(receipt_path, receipt)
    with pytest.raises(ValueError, match="Identical sequence windows"):
        ej.collect(path, "reported_baseline", True)
