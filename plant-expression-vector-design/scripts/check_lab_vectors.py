#!/usr/bin/env python3
"""Read-only, fresh per-invocation LAB_VECTOR_DIR probe. Standard library only."""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import socket
import stat
import sys
import uuid
from pathlib import Path

DEFAULT_CONFIG = Path(__file__).resolve().parents[1] / "templates/lab-library.config.json"
VARIABLE = "LAB_VECTOR_DIR"


def resolve_directory(config, override=None, environment=None):
    environment = os.environ if environment is None else environment
    if not isinstance(config, dict):
        raise ValueError("Library configuration must be a JSON object.")
    if config.get("variable_name") != VARIABLE:
        raise ValueError("Configuration must define LAB_VECTOR_DIR.")
    if config.get("read_only") is not True or config.get("check_every_invocation") is not True:
        raise ValueError("Source library must remain read-only and be checked every invocation.")
    if override is not None:
        value, origin = override, "explicit_argument"
    elif VARIABLE in environment:
        value, origin = environment[VARIABLE], "environment_variable"
    else:
        value, origin = config.get("default_path"), "skill_default"
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ValueError("LAB_VECTOR_DIR is empty or invalid; no silent fallback.")
    if not Path(value).is_absolute():
        raise ValueError("LAB_VECTOR_DIR must be an absolute directory path.")
    expected = config.get("expected_platform")
    if expected not in {"darwin", "linux", "win32", "any"}:
        raise ValueError("Declare the expected execution platform.")
    return value, origin, expected


def probe_directory(path, origin="skill_default", expected_platform="darwin", platform_name=None):
    observed = sys.platform if platform_name is None else platform_name
    result = {
        "schema_version": "1.0", "variable_name": VARIABLE, "path": str(path),
        "path_source": origin, "invocation_id": str(uuid.uuid4()),
        "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "execution_hostname": socket.gethostname(), "expected_platform": expected_platform,
        "observed_platform": observed, "status": "not_checked",
        "exists": None, "is_directory": None, "readable": None, "searchable": None,
        "directory_probe_verified": False, "inventory_performed": False,
        "source_read_only": True, "source_modified": False,
    }
    if expected_platform != "any" and expected_platform != observed:
        result.update(
            status="wrong_execution_host",
            next_action="Check on the authorized source device; this is not evidence the Mac directory is missing.",
        )
        return result
    try:
        info = os.stat(path)
        result["exists"] = True
        result["is_directory"] = stat.S_ISDIR(info.st_mode)
        if not result["is_directory"]:
            result.update(status="not_a_directory", next_action="Ask the user for a directory.")
            return result
        result["readable"] = os.access(path, os.R_OK)
        result["searchable"] = os.access(path, os.X_OK)
        if not result["readable"] or not result["searchable"]:
            result.update(
                status="permission_denied",
                next_action="Ask the user to confirm read/search access; do not change security settings.",
            )
            return result
        # One bounded directory access verifies the live mount without disclosing names.
        with os.scandir(path) as entries:
            next(entries, None)
        result.update(
            status="ready", directory_probe_verified=True,
            next_action="Use bounded, project-relevant inventory searches; do not modify the source library.",
        )
    except (FileNotFoundError, NotADirectoryError) as exc:
        result.update(
            status="missing_or_unavailable", exists=False, is_directory=False,
            directory_probe_verified=False, error_type=type(exc).__name__,
            next_action="Ask the user to confirm the path or mount; do not create a replacement directory.",
        )
    except PermissionError as exc:
        result.update(
            status="permission_denied", directory_probe_verified=False,
            error_type=type(exc).__name__,
            next_action="Access is denied; actual existence may be unconfirmed. Do not alter permissions.",
        )
    except OSError as exc:
        result.update(
            status="io_error", directory_probe_verified=False,
            error_type=type(exc).__name__, errno=exc.errno,
            next_action="Report the filesystem/connection failure without asserting the directory is absent.",
        )
    return result


def check(config_path=DEFAULT_CONFIG, override=None, environment=None):
    config = json.loads(Path(config_path).read_text(encoding="utf-8"))
    path, origin, expected = resolve_directory(config, override, environment)
    return probe_directory(path, origin, expected)


def save_report(path, report):
    path = Path(path)
    source = report.get("path")
    if not source:
        raise ValueError("Cannot write a report before the source path is resolved; use stdout diagnostics.")
    if source:
        try:
            path.resolve(strict=False).relative_to(Path(source).resolve(strict=False))
        except ValueError:
            pass
        else:
            raise ValueError("Do not write reports into LAB_VECTOR_DIR; the source library is read-only.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--lab-vector-dir", help="Explicit override for this invocation only.")
    parser.add_argument("--out", type=Path, help="New report file outside the source directory; never a cache.")
    args = parser.parse_args()
    try:
        report = check(args.config, args.lab_vector_dir)
    except (ValueError, TypeError, OSError) as exc:
        report = {
            "variable_name": VARIABLE, "status": "invalid_configuration",
            "checked_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "invocation_id": str(uuid.uuid4()), "error": str(exc),
            "exists": None, "source_modified": False,
        }
    if args.out:
        try:
            save_report(args.out, report)
        except (ValueError, OSError) as exc:
            print(json.dumps({"status": "report_not_written", "error": str(exc), "probe_result": report}, ensure_ascii=False))
            return 2
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "ready" else 2


if __name__ == "__main__":
    raise SystemExit(main())
