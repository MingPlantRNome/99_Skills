#!/usr/bin/env python3
"""Prepare explicit windows and run a verified, user-configured Evo2 command.

Standard library only. --weights-path belongs to THIS adapter, not EasyEvo2.
No SSH, scheduler submission, model download, shell evaluation or credentials.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import signal
import string
import subprocess
import time
from pathlib import Path

FIELDS = {"input", "output", "weights_path", "model_type", "window_size", "step_size", "device"}


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def read_fasta(path):
    records, name, seq = [], None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq).upper()))
            words = line[1:].split()
            if not words:
                raise ValueError("Empty FASTA ID.")
            name, seq = words[0], []
        else:
            if name is None:
                raise ValueError("FASTA sequence before header.")
            seq.append(line)
    if name is not None:
        records.append((name, "".join(seq).upper()))
    if not records or len({x[0] for x in records}) != len(records):
        raise ValueError("FASTA must contain unique, nonempty IDs.")
    for name, sequence in records:
        if len(sequence) < 2 or set(sequence) - set("ACGT"):
            raise ValueError(f"{name}: at least two unambiguous DNA bases required.")
    return records


def limit(profile, name):
    value = profile.get("limits", {}).get(name)
    if type(value) is not int or value <= 0:
        raise ValueError(f"Positive integer limits.{name} required.")
    return value


def fields_in(template):
    fields = set()
    for _, name, spec, conversion in string.Formatter().parse(template):
        if name is not None:
            if name not in FIELDS or spec or conversion:
                raise ValueError(f"Unsupported command placeholder: {name}")
            fields.add(name)
    return fields


def command(job, directory):
    profile = job["profile"]
    if profile.get("verified") is not True:
        raise ValueError("HPC command profile is unverified; inspect the installed loader first.")
    if profile.get("score_type") != "mean_log_likelihood":
        raise ValueError("Only explicitly verified mean_log_likelihood collection is supported.")
    if not profile.get("model_revision") or not profile.get("code_revision"):
        raise ValueError("Record model_revision and code_revision before execution.")
    context = profile.get("context_limit_nt")
    if type(context) is not int or context < max(x["length"] for x in job["items"]):
        raise ValueError("Actual model context_limit_nt must be verified and cover each input.")
    template = profile.get("argv_template")
    env_template = profile.get("env_template", {})
    if not isinstance(template, list) or not template or any(not isinstance(x, str) for x in template):
        raise ValueError("Provide a verified argv_template list; no guessed EasyEvo2 flags.")
    if not isinstance(env_template, dict) or any(not isinstance(v, str) for v in env_template.values()):
        raise ValueError("env_template must contain string values.")
    all_fields = set()
    for value in [*template, *env_template.values()]:
        all_fields |= fields_in(value)
    if not {"input", "output", "weights_path"} <= all_fields:
        raise ValueError("Verified command/environment must consume input, output and weights_path.")
    values = {
        "input": str((directory / "evo2_input.fasta").resolve()),
        "output": str((directory / "raw_scores.csv").resolve()),
        "weights_path": job["weights_path"], "model_type": profile["model_type"],
        "window_size": str(job["window_size"]), "step_size": str(job["step_size"]),
        "device": profile.get("device", "cuda:0"),
    }
    argv = [item.format_map(values) for item in template]
    if Path(argv[0]).name in {"sh", "bash", "zsh", "fish", "cmd", "powershell"}:
        raise ValueError("Shell interpreter templates are not supported; use a reviewed executable/script.")
    environment = {k: v.format_map(values) for k, v in env_template.items()}
    for key in environment:
        if any(s in key.upper() for s in ("TOKEN", "PASSWORD", "PRIVATE_KEY", "SECRET")):
            raise ValueError("Do not store credentials in an Evo2 command profile.")
    return argv, environment


def prepare(input_path, profile_path, weights_path, out, window_size=512, step_size=256):
    if not weights_path or not Path(weights_path).is_absolute() or Path(weights_path).suffix.lower() != ".pt":
        raise ValueError("Ask the user for the actual absolute HPC .pt path at runtime.")
    if type(window_size) is not int or window_size < 2 or not 1 <= step_size <= window_size:
        raise ValueError("Require window_size >= 2 and 1 <= step_size <= window_size.")
    profile = json.loads(Path(profile_path).read_text(encoding="utf-8"))
    if not isinstance(profile.get("model_type"), str) or not profile["model_type"]:
        raise ValueError("A model_type matching the actual checkpoint is required.")
    if window_size > limit(profile, "max_item_nt"):
        raise ValueError("Window size exceeds the configured resource limit.")
    records = read_fasta(input_path)
    items, lines, parents = [], [], []
    for parent_id, sequence in records:
        n = len(sequence)
        parents.append({"id": parent_id, "length": n, "sha256": hashlib.sha256(sequence.encode()).hexdigest()})
        if n <= window_size:
            starts = [0]
        else:
            starts = list(range(0, n - window_size + 1, step_size))
            if starts[-1] != n - window_size:
                starts.append(n - window_size)
        for start in starts:
            part = sequence[start:start + window_size]
            item_id = f"win{len(items)+1:06d}"
            lines += [f">{item_id}", part]
            items.append({
                "id": item_id, "parent_id": parent_id, "parent_length": n,
                "start": start + 1, "end": start + len(part), "length": len(part), "strand": 1,
                "sha256": hashlib.sha256(part.encode()).hexdigest(),
            })
    if len(items) > limit(profile, "max_items"):
        raise ValueError("Item count exceeds the approved resource limit.")
    limit(profile, "timeout_seconds")
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Job directory must be new or empty.")
    out.mkdir(parents=True, exist_ok=True)
    (out / "evo2_input.fasta").write_text("\n".join(lines) + "\n", encoding="utf-8")
    job = {
        "schema_version": "1.0", "created_at": now(),
        "weights_path": weights_path, "weights_path_scope": "HPC runtime; not checked on cloud",
        "profile": profile, "window_size": window_size, "step_size": step_size,
        "input_sha256": sha(out / "evo2_input.fasta"), "parents": parents, "items": items,
        "strand_policy": "provided_direction_only", "calibration_status": "not_calibrated",
        "execution_status": "not_executed",
    }
    dump(out / "job.json", job)
    return job


def load_job(path):
    path = Path(path).resolve()
    job = json.loads(path.read_text(encoding="utf-8"))
    directory = path.parent
    input_path = directory / "evo2_input.fasta"
    if sha(input_path) != job["input_sha256"]:
        raise ValueError("Input FASTA hash changed after preparation.")
    records = dict(read_fasta(input_path))
    if set(records) != {item["id"] for item in job["items"]} or len(records) != len(job["items"]):
        raise ValueError("Task IDs do not match the actual input FASTA.")
    for item in job["items"]:
        if hashlib.sha256(records[item["id"]].encode()).hexdigest() != item["sha256"]:
            raise ValueError("Per-item sequence hash mismatch.")
    return job, directory


def run(job_path, execute=False):
    job, directory = load_job(job_path)
    argv, environment = command(job, directory)
    if len(job["items"]) > limit(job["profile"], "max_items") or max(x["length"] for x in job["items"]) > limit(job["profile"], "max_item_nt"):
        raise ValueError("Task exceeds configured admission limits.")
    plan = {"status": "dry_run", "argv": argv, "env": environment, "weights_path": job["weights_path"]}
    if not execute:
        return plan
    weights = Path(job["weights_path"])
    if not weights.is_file():
        raise ValueError("The .pt file is not visible on this execution host; run on the authorized HPC.")
    output = directory / "raw_scores.csv"
    receipt_path = directory / "run_receipt.json"
    if output.exists() or receipt_path.exists():
        raise ValueError("Existing output/receipt: create a new versioned job directory.")
    stat = weights.stat()
    receipt = {
        "status": "running", "started_at": now(), "input_sha256": job["input_sha256"],
        "job_sha256": sha(job_path), "argv": argv, "profile_id": job["profile"].get("profile_id"),
        "weights_fingerprint": {"method": "path_size_mtime_not_full_hash", "path": str(weights), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns},
        "simulation": bool(job["profile"].get("simulation", False)),
        "biological_validation": "not_performed",
    }
    dump(receipt_path, receipt)
    env = os.environ.copy()
    env.update(environment)
    started = time.monotonic()
    proc = None
    try:
        with (directory / "run.log").open("wb") as handle:
            proc = subprocess.Popen(
                argv, cwd=directory, env=env, stdout=handle, stderr=subprocess.STDOUT,
                shell=False, start_new_session=True,
            )
            try:
                code = proc.wait(timeout=limit(job["profile"], "timeout_seconds"))
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                raise TimeoutError("Configured job timeout reached; inspect scheduler state if a launcher was used.")
        receipt["returncode"] = code
        if code != 0 or not output.is_file():
            raise RuntimeError(f"Command failed or output missing (returncode={code}).")
        receipt["output_sha256"] = sha(output)
        receipt["status"] = "simulation_completed" if receipt["simulation"] else "completed"
    except Exception as exc:
        receipt["status"] = "failed"
        receipt["error"] = str(exc)
    finally:
        receipt["finished_at"] = now()
        receipt["elapsed_seconds"] = time.monotonic() - started
        dump(receipt_path, receipt)
    if receipt["status"] == "failed":
        raise RuntimeError(receipt["error"])
    return receipt


def collect(job_path, baseline_id, allow_simulation=False):
    job, directory = load_job(job_path)
    command(job, directory)
    receipt_path = directory / "run_receipt.json"
    if not receipt_path.exists():
        raise ValueError("No actual run receipt; do not fabricate model execution evidence.")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("simulation") and not allow_simulation:
        raise ValueError("Simulation is not Evo2 inference; test-only opt-in required.")
    if receipt.get("status") not in {"completed", "simulation_completed"}:
        raise ValueError("Actual execution did not complete.")
    if receipt.get("input_sha256") != job["input_sha256"] or receipt.get("job_sha256") != sha(job_path):
        raise ValueError("Run receipt does not match the job/input.")
    raw = directory / "raw_scores.csv"
    if sha(raw) != receipt.get("output_sha256"):
        raise ValueError("Raw output changed after execution.")
    profile = job["profile"]
    columns = profile["score_columns"]
    delimiter = profile.get("output_delimiter", ",")
    if not isinstance(delimiter, str) or len(delimiter) != 1:
        raise ValueError("Output delimiter must be a single character.")
    wanted = {x["id"] for x in job["items"]}
    scores, extra_rows = {}, 0
    with raw.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not {columns["id"], columns["value"]} <= set(reader.fieldnames or []):
            raise ValueError("Actual output is missing configured score columns.")
        for row in reader:
            name = row[columns["id"]]
            if name not in wanted:
                extra_rows += 1
                continue
            value = float(row[columns["value"]])
            if name in scores or not math.isfinite(value):
                raise ValueError("Duplicate or nonfinite model score.")
            scores[name] = value
    if set(scores) != wanted:
        raise ValueError("Missing expected model scores.")
    tolerance = profile.get("score_tolerance", 1e-5)
    if not isinstance(tolerance, (float, int)) or not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("score_tolerance must be a finite nonnegative number.")
    by_hash = {}
    for item in job["items"]:
        value = scores[item["id"]]
        previous = by_hash.setdefault(item["sha256"], value)
        if not math.isclose(value, previous, rel_tol=0, abs_tol=tolerance):
            raise ValueError("Identical sequence windows have inconsistent scores; investigate batch/padding/protocol.")
    baseline = [x for x in job["items"] if x["parent_id"] == baseline_id]
    if not baseline:
        raise ValueError("baseline-id is not a prepared parent sequence.")
    references = {(x["parent_length"], x["start"], x["end"], x["strand"]): scores[x["id"]] for x in baseline}
    rows = []
    for item in job["items"]:
        key = (item["parent_length"], item["start"], item["end"], item["strand"])
        ref = references.get(key)
        rows.append({
            **item, "model_type": profile["model_type"], "score_type": profile["score_type"],
            "score": scores[item["id"]], "baseline_score": ref,
            "delta_score": scores[item["id"]] - ref if ref is not None else "NA",
            "comparability": "paired_window_only" if ref is not None else "incomparable_length_or_coordinates",
            "execution_status": receipt["status"], "calibration_status": "not_calibrated",
        })
    out = directory / "evo2_scores.tsv"
    if out.exists():
        raise ValueError("Collected output already exists; preserve previous results.")
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "execution_status": receipt["status"], "score_rows": len(rows), "ignored_extra_rows": extra_rows,
        "input_sha256": job["input_sha256"], "model_type": profile["model_type"],
        "score_type": profile["score_type"], "calibration_status": "not_calibrated",
        "interpretation": "Model preference only; no inferred plant expression, splicing or stability improvement.",
    }
    dump(directory / "collection_report.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--profile", type=Path, required=True)
    p.add_argument("--weights-path", required=True, help="Runtime absolute .pt path supplied by the user; not an EasyEvo2 flag.")
    p.add_argument("--window-size", type=int, default=512)
    p.add_argument("--step-size", type=int, default=256)
    p.add_argument("--out", type=Path, required=True)
    p = commands.add_parser("run")
    p.add_argument("--job", type=Path, required=True)
    p.add_argument("--execute", action="store_true")
    p = commands.add_parser("collect")
    p.add_argument("--job", type=Path, required=True)
    p.add_argument("--baseline-id", required=True)
    p.add_argument("--allow-simulation", action="store_true", help="Unit tests only; never represents Evo2 inference.")
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            result = prepare(args.input, args.profile, args.weights_path, args.out, args.window_size, args.step_size)
            print(json.dumps({"status": "prepared_not_executed", "items": len(result["items"]), "job": str(args.out / "job.json")}))
        elif args.command == "run":
            print(json.dumps(run(args.job, args.execute), ensure_ascii=False))
        else:
            print(json.dumps(collect(args.job, args.baseline_id, args.allow_simulation), ensure_ascii=False))
        return 0
    except Exception as exc:
        print(json.dumps({"status": "fail", "error": str(exc)}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
