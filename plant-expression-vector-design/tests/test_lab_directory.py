"""Temporary-directory and simulated-platform tests; no real NAS inventory."""
from __future__ import annotations

import copy
import errno
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("check_lab_vectors", ROOT / "scripts/check_lab_vectors.py")
lb = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lb)


def config(path):
    return {
        "variable_name": "LAB_VECTOR_DIR", "default_path": str(path),
        "expected_platform": "darwin", "read_only": True, "check_every_invocation": True,
    }


def test_default_configuration_is_single_source():
    data = json.loads((ROOT / "templates/lab-library.config.json").read_text())
    path, origin, platform = lb.resolve_directory(data, environment={})
    assert origin == "skill_default" and platform == "darwin"
    assert path.endswith("/403_Plasmid/Maps")
    assert data["check_every_invocation"] and data["read_only"]


def test_explicit_environment_default_priority(tmp_path):
    data = config(tmp_path / "default")
    env = {"LAB_VECTOR_DIR": str(tmp_path / "environment")}
    assert lb.resolve_directory(data, str(tmp_path / "explicit"), env)[:2] == (
        str(tmp_path / "explicit"), "explicit_argument"
    )
    assert lb.resolve_directory(data, environment=env)[:2] == (
        str(tmp_path / "environment"), "environment_variable"
    )
    assert lb.resolve_directory(data, environment={})[:2] == (
        str(tmp_path / "default"), "skill_default"
    )


@pytest.mark.parametrize("value", ["", "   ", "relative/path", "\x00invalid"])
def test_invalid_explicit_override_never_falls_back(tmp_path, value):
    with pytest.raises(ValueError):
        lb.resolve_directory(config(tmp_path), value, {"LAB_VECTOR_DIR": str(tmp_path)})


def test_empty_environment_never_silently_uses_default(tmp_path):
    with pytest.raises(ValueError, match="empty"):
        lb.resolve_directory(config(tmp_path), environment={"LAB_VECTOR_DIR": ""})


def test_configuration_cannot_disable_preflight_or_readonly(tmp_path):
    for field in ("read_only", "check_every_invocation"):
        data = config(tmp_path)
        data[field] = False
        with pytest.raises(ValueError, match="read-only"):
            lb.resolve_directory(data, environment={})
    with pytest.raises(ValueError, match="JSON object"):
        lb.resolve_directory([], environment={})


def test_ready_does_not_inventory_or_write_source(tmp_path):
    source = tmp_path / "library with spaces"
    source.mkdir()
    file = source / "do-not-read.fa"
    file.write_text("Private test placeholder, not a sequence.")
    before = {p.name: (p.stat().st_mtime_ns, p.read_bytes()) for p in source.iterdir()}
    result = lb.probe_directory(source, platform_name="darwin")
    after = {p.name: (p.stat().st_mtime_ns, p.read_bytes()) for p in source.iterdir()}
    assert result["status"] == "ready"
    assert result["exists"] is True and result["directory_probe_verified"]
    assert result["inventory_performed"] is False and result["source_modified"] is False
    assert before == after
    assert "do-not-read" not in json.dumps(result)


def test_empty_directory_is_ready_not_plasmid_found(tmp_path):
    result = lb.probe_directory(tmp_path, platform_name="darwin")
    assert result["status"] == "ready" and result["inventory_performed"] is False


def test_missing_directory_is_not_created(tmp_path):
    source = tmp_path / "offline-or-missing"
    result = lb.probe_directory(source, platform_name="darwin")
    assert result["status"] == "missing_or_unavailable" and result["exists"] is False
    assert not source.exists()


def test_file_is_not_a_directory(tmp_path):
    file = tmp_path / "a-file"
    file.write_text("fixture")
    result = lb.probe_directory(file, platform_name="darwin")
    assert result["status"] == "not_a_directory" and result["exists"] is True


def test_wrong_host_does_not_probe_or_claim_missing(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Must not touch a Mac path on wrong host.")
    monkeypatch.setattr(lb.os, "stat", forbidden)
    result = lb.probe_directory("/Users/fixture/Maps", platform_name="linux")
    assert result["status"] == "wrong_execution_host" and result["exists"] is None


def test_permission_failure_does_not_assert_missing(tmp_path, monkeypatch):
    def denied(*args, **kwargs):
        raise PermissionError("simulated permission failure")
    monkeypatch.setattr(lb.os, "stat", denied)
    result = lb.probe_directory(tmp_path, platform_name="darwin")
    assert result["status"] == "permission_denied" and result["exists"] is None


def test_read_access_denied(tmp_path, monkeypatch):
    monkeypatch.setattr(lb.os, "access", lambda *args: False)
    result = lb.probe_directory(tmp_path, platform_name="darwin")
    assert result["status"] == "permission_denied" and result["exists"] is True


def test_io_failure_is_not_missing(tmp_path, monkeypatch):
    def failed(*args, **kwargs):
        raise OSError(errno.EIO, "simulated filesystem failure")
    monkeypatch.setattr(lb.os, "stat", failed)
    result = lb.probe_directory(tmp_path, platform_name="darwin")
    assert result["status"] == "io_error" and result["exists"] is None


def test_scandir_permission_race_reports_denied(tmp_path, monkeypatch):
    def denied(*args, **kwargs):
        raise PermissionError("simulated access changed after stat")
    monkeypatch.setattr(lb.os, "scandir", denied)
    result = lb.probe_directory(tmp_path, platform_name="darwin")
    assert result["status"] == "permission_denied" and not result["directory_probe_verified"]


def test_each_invocation_rechecks_instead_of_caching(tmp_path):
    source = tmp_path / "transient"
    source.mkdir()
    first = lb.probe_directory(source, platform_name="darwin")
    source.rmdir()
    second = lb.probe_directory(source, platform_name="darwin")
    assert first["status"] == "ready" and second["status"] == "missing_or_unavailable"
    assert first["invocation_id"] != second["invocation_id"]


def test_report_cannot_be_saved_in_source(tmp_path):
    report = lb.probe_directory(tmp_path, platform_name="darwin")
    with pytest.raises(ValueError, match="read-only"):
        lb.save_report(tmp_path / "new-folder/report.json", report)
    assert not (tmp_path / "new-folder").exists()


def test_report_is_new_and_not_overwritten(tmp_path):
    source = tmp_path / "library"
    source.mkdir()
    report = lb.probe_directory(source, platform_name="darwin")
    output = tmp_path / "project/preflight.json"
    lb.save_report(output, report)
    assert json.loads(output.read_text())["invocation_id"] == report["invocation_id"]
    with pytest.raises(FileExistsError):
        lb.save_report(output, report)


def test_cli_default_on_cloud_is_not_missing_directory():
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts/check_lab_vectors.py")],
        capture_output=True, text=True, env={},
    )
    if sys.platform != "darwin":
        assert result.returncode == 2
        payload = json.loads(result.stdout)
        assert payload["status"] == "wrong_execution_host" and payload["exists"] is None


def test_symlink_report_escape_is_blocked(tmp_path):
    source = tmp_path / "library"
    source.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(source, target_is_directory=True)
    report = lb.probe_directory(source, platform_name="darwin")
    with pytest.raises(ValueError, match="read-only"):
        lb.save_report(alias / "report.json", report)
    assert list(source.iterdir()) == []


def test_unresolved_source_cannot_bypass_readonly_report_guard(tmp_path):
    target = tmp_path / "unknown-source/report.json"
    with pytest.raises(ValueError, match="source path is resolved"):
        lb.save_report(target, {"status": "invalid_configuration"})
    assert not target.parent.exists()
