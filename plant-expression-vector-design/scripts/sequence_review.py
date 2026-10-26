#!/usr/bin/env python3
"""Read-only observations on annotated 5' leaders and exact regulatory repeats.

Not an expression, splicing, silencing or recombination predictor. No sequence edits.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition


def role_of(feature):
    if feature.type in {"5'UTR", "5UTR", "five_prime_UTR"}:
        return "leader"
    if feature.type in {"promoter", "terminator"}:
        return feature.type
    for key in ("regulatory_class", "design_role"):
        value = feature.qualifiers.get(key, [])
        if value and value[0] in {"promoter", "terminator", "leader"}:
            return value[0]
    return None


def identity(feature, index):
    value = feature.qualifiers.get("design_id", feature.qualifiers.get("label", []))
    return value[0] if value else f"annotation_{index}"


def review(record, min_repeat_nt=20):
    sequence = str(record.seq).upper()
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("Review requires an exact nonempty DNA record; ambiguous bases are not silently ignored.")
    if type(min_repeat_nt) is not int or min_repeat_nt < 1:
        raise ValueError("min_repeat_nt is an algorithmic size filter, not a biological threshold.")
    leaders, repeats, flags = [], {}, []
    observed_loci, observed_leaders = set(), set()
    for index, feature in enumerate(record.features, 1):
        role = role_of(feature)
        if role is None:
            continue
        fid = identity(feature, index)
        if feature.location is None:
            flags.append({"category": "unassessed_annotation", "feature": fid, "observation": "Missing feature location."})
            continue
        for part in feature.location.parts:
            if not isinstance(part.start, ExactPosition) or not isinstance(part.end, ExactPosition):
                raise ValueError(f"{fid}: resolve fuzzy boundaries before coordinate-based review.")
            if not 0 <= int(part.start) < int(part.end) <= len(record):
                raise ValueError(f"{fid}: feature boundary out of range.")
        if feature.location.strand not in {1, -1}:
            flags.append({"category": "unassessed_annotation", "feature": fid, "observation": "Direction unknown; no strand was guessed."})
            continue
        extracted = str(feature.extract(Seq(sequence))).upper()
        positions = list(feature.location)
        # Same physical annotation repeated under multiple labels is not a second DNA copy.
        locus = tuple(sorted(positions))
        is_new_regulatory_locus = locus not in observed_loci
        observed_loci.add(locus)
        leader_key = (locus, feature.location.strand)
        if role == "leader" and leader_key not in observed_leaders:
            observed_leaders.add(leader_key)
            candidates = []
            for start in range(len(extracted) - 2):
                if extracted[start:start + 3] != "ATG":
                    continue
                stop = next((
                    pos for pos in range(start + 3, len(extracted) - 2, 3)
                    if extracted[pos:pos + 3] in {"TAA", "TAG", "TGA"}
                ), None)
                candidate = {
                    "leader_start_nt": start + 1,
                    "start_triplet_positions_1based": [p + 1 for p in positions[start:start + 3]],
                    "context": extracted[max(0, start - 3):min(len(extracted), start + 6)],
                    "candidate_class": "ATG_ORF_within_annotated_leader" if stop is not None else "ATG_candidate_stop_not_resolved_within_leader",
                    "stop_leader_start_nt": stop + 1 if stop is not None else None,
                    "length_aa_before_stop": (stop - start) // 3 if stop is not None else None,
                    "translation_or_function_validated": False,
                }
                candidates.append(candidate)
            leaders.append({"feature": fid, "strand": feature.location.strand, "length_nt": len(extracted), "candidates": candidates})
            if candidates:
                flags.append({
                    "category": "leader_ATG_candidates", "feature": fid,
                    "observation": f"{len(candidates)} exact ATG candidates in the annotated direction; not a functional uORF call.",
                })
        if is_new_regulatory_locus and len(extracted) >= min_repeat_nt:
            reverse = str(Seq(extracted).reverse_complement())
            canonical = min(extracted, reverse)
            key = hashlib.sha256(canonical.encode()).hexdigest()
            repeats.setdefault(key, []).append({
                "feature": fid, "role": role, "length_nt": len(extracted),
                "segments_1based": [[int(p.start) + 1, int(p.end)] for p in feature.location.parts],
                "strand": feature.location.strand,
            })
    repeated = [{"canonical_sequence_sha256": key, "copies": copies} for key, copies in repeats.items() if len(copies) > 1]
    for group in repeated:
        flags.append({
            "category": "exact_annotated_regulatory_repeat",
            "feature": "; ".join(x["feature"] for x in group["copies"]),
            "observation": f"{len(group['copies'])} distinct annotated loci have the same exact sequence up to reverse complement; no silencing/recombination outcome inferred.",
        })
    if not leaders:
        flags.append({"category": "not_assessed", "feature": "5prime_leader", "observation": "No readable directional 5'UTR annotation; no leader was inferred from a promoter name."})
    if not repeats:
        flags.append({"category": "not_assessed", "feature": "regulatory_repeats", "observation": "No eligible annotated regulatory sequence above the stated size filter."})
    return {
        "record_id": record.id, "length_nt": len(record),
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "status": "sequence_observations_only", "source_modified": False,
        "leader_review": leaders, "exact_regulatory_repeats": repeated,
        "flags": flags, "min_repeat_nt": min_repeat_nt,
        "limitations": [
            "ATG-only scan; non-AUG initiation, RNA structure and actual translation are not evaluated.",
            "ORFs are only resolved within the supplied leader annotation; overlap into the main CDS is not inferred.",
            "Only exact complete annotated regulatory sequences up to reverse complement are compared; no approximate/global homology scan.",
            "No expression, RNA-processing, stability, silencing, recombination or risk probability prediction.",
            "Observations do not authorize deleting starts, recoding sequences or replacing reported parts.",
        ],
    }


def generate(input_path, out, record_id=None, min_repeat_nt=20):
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new or empty versioned review output directory.")
    records = list(SeqIO.parse(input_path, "genbank"))
    if record_id is not None:
        records = [r for r in records if r.id == record_id]
    if len(records) != 1:
        raise ValueError("Select one exact GenBank record using --record-id if necessary.")
    result = review(records[0], min_repeat_nt)
    out.mkdir(parents=True, exist_ok=True)
    (out / "sequence_review.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (out / "review_flags.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["category", "feature", "observation"], delimiter="\t")
        writer.writeheader()
        writer.writerows(result["flags"])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--record-id")
    parser.add_argument("--min-repeat-nt", type=int, default=20)
    args = parser.parse_args()
    try:
        result = generate(args.input, args.out, args.record_id, args.min_repeat_nt)
        print(json.dumps({"status": result["status"], "flags": len(result["flags"]), "source_modified": False}))
        return 0
    except Exception as exc:
        print(json.dumps({"status": "fail", "error": str(exc)}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
