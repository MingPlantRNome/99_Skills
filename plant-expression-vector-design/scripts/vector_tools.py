#!/usr/bin/env python3
"""Offline, evidence-gated DNA assembly/export and equal-length candidate checks."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import html
import json
import tempfile
from pathlib import Path

import jsonschema
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, ExactPosition, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[1]


def dump_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def digest(sequence):
    return hashlib.sha256(str(sequence).upper().encode("ascii")).hexdigest()


def dna(value):
    value = str(value).upper()
    if not value or set(value) - set("ACGT"):
        raise ValueError("Exact DNA must be nonempty ACGT; no N, RNA U or placeholders.")
    return value


def manifest(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    schema = json.loads((ROOT / "templates/construct.schema.json").read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(data)
    for key in ("sources", "features", "protected_regions"):
        ids = [x["id"] for x in data.get(key, [])]
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate {key} IDs.")
    return data


def read_record(path, fmt=None, record_id=None):
    path = Path(path)
    if fmt is None:
        fmt = "genbank" if path.suffix.lower() in {".gb", ".gbk", ".genbank"} else "fasta"
    records = list(SeqIO.parse(path, fmt))
    if record_id is not None:
        records = [r for r in records if r.id == record_id]
    if len(records) != 1:
        raise ValueError(f"{path.name}: expected one exact record; specify a unique record_id.")
    record = records[0]
    record.seq = Seq(dna(record.seq))
    return record


def location(spec, length):
    parts = []
    seen = set()
    for start, end in spec["segments"]:
        if not (1 <= start <= end <= length):
            raise ValueError(f"{spec['id']}: coordinates outside 1..{length}; split origin-spanning segments.")
        positions = set(range(start, end + 1))
        if seen & positions:
            raise ValueError(f"{spec['id']}: overlapping joined segments.")
        seen.update(positions)
        parts.append(SimpleLocation(start - 1, end, strand=spec["strand"]))
    return parts[0] if len(parts) == 1 else CompoundLocation(parts, operator="join")


def explicit_feature(spec, length):
    quals = {
        k: [v] if isinstance(v, str) else list(v)
        for k, v in spec.get("qualifiers", {}).items()
    }
    quals["design_id"] = [spec["id"]]
    if spec["type"] == "CDS":
        quals["translation"] = [spec["expected_protein"].rstrip("*")]
        quals["design_complete_cds"] = ["true" if spec["complete_cds"] else "false"]
        quals["transl_table"] = [str(spec.get("translation_table", 1))]
    return SeqFeature(location(spec, length), type=spec["type"], id=spec["id"], qualifiers=quals)


def feature_id(feature):
    return feature.qualifiers.get("design_id", [feature.id])[0]


def protein(feature, sequence):
    if "codon_start" in feature.qualifiers and feature.qualifiers["codon_start"] != ["1"]:
        raise ValueError(f"{feature_id(feature)}: partial codon_start requires specialist reannotation.")
    complete = feature.qualifiers.get("design_complete_cds", ["true"])[0] == "true"
    table = int(feature.qualifiers.get("transl_table", ["1"])[0])
    coding = feature.extract(Seq(sequence))
    if len(coding) % 3:
        raise ValueError(f"{feature_id(feature)}: CDS length not divisible by 3.")
    translated = str(coding.translate(table=table, cds=complete))
    if "*" in translated.rstrip("*"):
        raise ValueError(f"{feature_id(feature)}: internal stop codon.")
    return translated.rstrip("*")


def check_features(record):
    results = []
    for f in record.features:
        for part in f.location.parts:
            if not isinstance(part.start, ExactPosition) or not isinstance(part.end, ExactPosition):
                raise ValueError(f"{feature_id(f)}: fuzzy annotation needs review.")
            if not (0 <= int(part.start) < int(part.end) <= len(record)):
                raise ValueError(f"{feature_id(f)}: invalid feature boundary.")
        if f.type == "CDS":
            if f.location.strand not in {1, -1}:
                raise ValueError(f"{feature_id(f)}: unresolved CDS direction.")
            expected = f.qualifiers.get("translation")
            if not expected:
                raise ValueError(f"{feature_id(f)}: no expected protein; do not infer an unverified CDS.")
            actual = protein(f, record.seq)
            if actual != expected[0].replace(" ", "").replace("\n", "").rstrip("*"):
                raise ValueError(f"{feature_id(f)}: expected protein mismatch.")
            results.append({"feature": feature_id(f), "translation": actual, "status": "pass"})
    return results


def site_hits(sequence, motif, circular=False):
    sequence, motif = dna(sequence), dna(motif)
    n = len(sequence)
    strands = [(motif, 1)]
    rc = str(Seq(motif).reverse_complement())
    if rc != motif:
        strands.append((rc, -1))
    if circular:
        seq = (sequence * ((n + len(motif) - 2) // n + 1))[:n + len(motif) - 1]
    else:
        seq = sequence
    hits = []
    for pattern, strand in strands:
        start = 0
        while True:
            found = seq.find(pattern, start)
            if found < 0 or found >= n:
                break
            hits.append({"start": found + 1, "strand": strand, "crosses_origin": found + len(motif) > n})
            start = found + 1
    return sorted(hits, key=lambda x: (x["start"], x["strand"]))


def sites_report(data, sequence):
    reports = []
    for site in data.get("restriction_sites", []):
        hits = site_hits(sequence, site["sequence"], data["topology"] == "circular")
        if hits and site["must_be_absent"]:
            raise ValueError(f"Forbidden site {site['name']} remains ({len(hits)} hits).")
        reports.append({**site, "hits": hits})
    return reports


def equivalent(candidate, reference, circular):
    candidate, reference = dna(candidate), dna(reference)
    if len(candidate) != len(reference):
        return False
    if not circular:
        return candidate == reference
    return candidate in reference + reference or str(Seq(candidate).reverse_complement()) in reference + reference


def fusion_source(spec, source_records, label):
    source_id = spec["source_id"]
    if source_id not in source_records:
        raise ValueError(f"{label}: unknown fusion source {source_id}.")
    source = source_records[source_id]
    loc = location({
        "id": label, "segments": spec["source_segments"], "strand": spec["source_strand"]
    }, len(source))
    return dna(loc.extract(source.seq))


def check_fusion_junctions(decision, target, record, source_records):
    """Check continuous ORF composition and source stops; no functional prediction."""
    plan = decision.get("fusion_plan")
    if not isinstance(plan, dict) or plan.get("assembly_mode") != "continuous_orf":
        raise ValueError("New tag fusion requires a complete continuous_orf fusion_plan.")
    if target.qualifiers.get("design_complete_cds", ["true"])[0] != "true":
        raise ValueError("New tag fusion requires a complete fused CDS.")
    table = int(target.qualifiers.get("transl_table", ["1"])[0])
    stop_codons = set(CodonTable.unambiguous_dna_by_id[table].stop_codons)
    coding = dna(target.extract(record.seq))
    if len(coding) % 3 or coding[-3:] not in stop_codons:
        raise ValueError("Complete fused CDS requires an in-frame final stop codon.")
    final = plan["final_stop"]
    if final["codon"] != coding[-3:] or fusion_source(final, source_records, "final_stop") != coding[-3:]:
        raise ValueError("Final stop codon does not match the declared source and fused CDS.")
    if not (plan.get("linker_flexibility_basis") or "").strip() or not plan.get("linker_evidence_ids"):
        raise ValueError("Flexible linker requires a declared rationale and evidence IDs.")
    positions = list(target.location)
    if len(positions) != len(set(positions)):
        raise ValueError("Fusion CDS must use nonoverlapping ordered coding segments.")
    offsets = {pos: index for index, pos in enumerate(positions)}
    all_features = {feature_id(f): f for f in record.features}
    components = plan["components"]
    ids = [c["feature_id"] for c in components]
    if len(ids) != len(set(ids)):
        raise ValueError("Fusion components must have separate unique feature IDs.")
    cursor, audit, tag_ids, major, junctions = 0, [], set(), [], []
    for index, component in enumerate(components):
        component_id = component["feature_id"]
        feature = all_features.get(component_id)
        if feature is None or feature is target:
            raise ValueError(f"{component_id}: missing separate fusion component annotation.")
        if feature.location.strand != target.location.strand:
            raise ValueError(f"{component_id}: fusion component is on the wrong strand.")
        mapped = [offsets.get(pos, -1) for pos in feature.location]
        if not mapped or mapped != list(range(cursor, cursor + len(mapped))) or len(mapped) % 3 or cursor % 3:
            raise ValueError(f"{component_id}: fusion component order, coverage or codon boundary mismatch.")
        actual = dna(feature.extract(record.seq))
        source = fusion_source(component, source_records, component_id)
        if len(source) % 3:
            raise ValueError(f"{component_id}: source coding segment is not in frame.")
        source_stop = source[-3:] if source[-3:] in stop_codons else None
        action = component["stop_action"]
        if source_stop is None:
            if action != "none_present":
                raise ValueError(f"{component_id}: no source terminal stop exists; cannot claim removal/retention.")
            body = source
        else:
            if action == "remove_for_fusion":
                body = source[:-3]
            elif action == "retained_as_final" and index == len(components) - 1:
                body = source[:-3]
                if final["provider"] != "last_component" or final["codon"] != source_stop:
                    raise ValueError(f"{component_id}: retained final stop conflicts with its provider/codon.")
            else:
                raise ValueError(f"{component_id}: upstream terminal stop must be removed for continuous fusion.")
        if actual != body:
            raise ValueError(f"{component_id}: final component DNA differs from its source after stop handling.")
        peptide = str(Seq(actual).translate(table=table))
        if "*" in peptide:
            raise ValueError(f"{component_id}: internal stop codon in fusion component.")
        if cursor == 0 and actual[:3] in CodonTable.unambiguous_dna_by_id[table].start_codons:
            peptide = "M" + peptide[1:]
        if peptide != component["expected_peptide"]:
            raise ValueError(f"{component_id}: component/linker peptide does not match the declared sequence.")
        role = component["role"]
        if role in {"target", "tag"}:
            major.append((index, role))
        if role == "tag":
            tag_ids.add(component_id)
        feature.qualifiers["design_role"] = [role]
        audit.append({
            "feature_id": component_id, "role": role, "source_terminal_stop": source_stop,
            "stop_action": action, "encoded_peptide": peptide, "length_aa": len(peptide),
            "coding_start": cursor + 1, "coding_end": cursor + len(actual),
            "linker_dna": actual if role == "flexible_linker" else None,
            "source_id": component["source_id"], "source_segments": component["source_segments"],
            "source_strand": component["source_strand"], "status": "pass",
        })
        cursor += len(actual)
        if index < len(components) - 1:
            start, end = max(0, cursor - 9), min(len(coding) - 3, cursor + 9)
            context_peptide = str(Seq(coding[start:end]).translate(table=table))
            if start == 0 and coding[:3] in CodonTable.unambiguous_dna_by_id[table].start_codons:
                context_peptide = "M" + context_peptide[1:]
            junctions.append({
                "upstream_feature_id": component_id,
                "downstream_feature_id": components[index + 1]["feature_id"],
                "coding_boundary_after_nt": cursor, "coding_start": start + 1, "coding_end": end,
                "dna_context": coding[start:end], "peptide_context": context_peptide,
                "coordinate_system": "1-based inclusive within spliced coding orientation, not genomic coordinates",
            })
    if cursor != len(coding) - 3:
        raise ValueError("Fusion components must cover the entire coding body; no hidden scar, overlap or gap.")
    if tag_ids != set(decision["tag_feature_ids"]) or not any(role == "target" for _, role in major):
        raise ValueError("Fusion target/tag components do not match the declared tag plan.")
    if final["provider"] == "last_component" and components[-1]["stop_action"] != "retained_as_final":
        raise ValueError("Final stop claimed from last component but not recorded as retained_as_final.")
    for (left, role_left), (right, role_right) in zip(major, major[1:]):
        if role_left != role_right and not any(
            c["role"] == "flexible_linker" for c in components[left + 1:right]
        ):
            raise ValueError("A flexible linker is required at every target/tag junction.")
    compressed = []
    for _, role in major:
        if not compressed or compressed[-1] != role:
            compressed.append(role)
    site = decision["fusion_site"]
    valid_site = (
        site == "C_terminal" and compressed == ["target", "tag"]
        or site == "N_terminal" and compressed == ["tag", "target"]
        or site == "internal" and len(compressed) >= 3
        and compressed[0] == compressed[-1] == "target"
    )
    if not valid_site:
        raise ValueError("Declared fusion site disagrees with the N-to-C component order.")
    return {
        "status": "continuous_orf_and_junctions_pass", "components": audit, "junctions": junctions,
        "final_stop": final, "flexibility_and_function": "not_experimentally_validated",
    }


def check_tag_decisions(data, record, source_records=None):
    """Validate declared choices/annotations, not user consent or literature truth."""
    decisions = data.get("tag_decisions", [])
    de_novo = data["replication_class"] == "de_novo"
    required_targets = set(data.get("new_expression_target_ids", []))
    if de_novo and not required_targets:
        raise ValueError("de_novo requires declared new expression CDS targets and their tag decisions.")
    by_target = {}
    features = {feature_id(f): f for f in record.features}
    for target_id in required_targets:
        if target_id not in features or features[target_id].type != "CDS":
            raise ValueError(f"{target_id}: declared new expression target must be an actual CDS.")
    checks = []
    source_records = {} if source_records is None else source_records
    for decision in decisions:
        target_id = decision["target_feature_id"]
        if target_id in by_target:
            raise ValueError(f"{target_id}: duplicate tag decisions.")
        by_target[target_id] = decision
        target = features.get(target_id)
        if target is None or target.type != "CDS":
            raise ValueError(f"{target_id}: tag decision target must be an actual CDS.")
        new_design = de_novo or target_id in required_targets or decision["design_context"] == "new_design"
        if decision["choice"] == "pending":
            raise ValueError(f"{target_id}: tag decision pending; ask the user before final export.")
        if new_design and (
            decision["decision_source"] != "user_explicit"
            or not (decision.get("user_response") or "").strip()
        ):
            raise ValueError(f"{target_id}: new design requires the actual explicit user tag choice.")
        if decision["design_context"] == "replicated_original" and decision["decision_source"] not in {
            "reported_design", "user_explicit"
        }:
            raise ValueError(f"{target_id}: original tag plan requires traceable source evidence.")
        tag_ids = decision["tag_feature_ids"]
        fusion_check = "not_applicable"
        if len(tag_ids) != len(set(tag_ids)):
            raise ValueError(f"{target_id}: duplicate tag feature IDs.")
        if decision["choice"] == "no_tag":
            if decision["selected_tag"] is not None or decision["fusion_site"] != "none" or tag_ids:
                raise ValueError(f"{target_id}: no_tag conflicts with declared tag/fusion annotations.")
        elif decision["choice"] == "add_tag":
            if decision["plan_confirmed"] is not True:
                raise ValueError(f"{target_id}: selected tag plan is not confirmed.")
            for field in ("purpose", "selected_tag", "sequence_source"):
                if not (decision.get(field) or "").strip():
                    raise ValueError(f"{target_id}: add_tag requires {field}.")
            if not decision["literature_evidence_ids"]:
                raise ValueError(f"{target_id}: add_tag requires literature evidence IDs, not a generic tag list.")
            if decision["fusion_site"] not in {"N_terminal", "C_terminal", "internal"} or not tag_ids:
                raise ValueError(f"{target_id}: add_tag requires a selected fusion site and tag feature IDs.")
            target_positions = {
                p for part in target.location.parts for p in range(int(part.start), int(part.end))
            }
            for tag_id in tag_ids:
                tag = features.get(tag_id)
                if tag is None or tag_id == target_id:
                    raise ValueError(f"{target_id}: missing or invalid separate tag feature {tag_id}.")
                tag_positions = {
                    p for part in tag.location.parts for p in range(int(part.start), int(part.end))
                }
                if not tag_positions <= target_positions or tag.location.strand != target.location.strand:
                    raise ValueError(f"{target_id}: tag feature must lie within the CDS coding segments on the same strand.")
            if new_design or decision.get("fusion_plan") is not None:
                fusion_check = check_fusion_junctions(decision, target, record, source_records)
            else:
                fusion_check = "not_audited_legacy_replication"
        else:
            raise ValueError(f"{target_id}: unsupported tag choice.")
        checks.append({
            "target_feature_id": target_id, "choice": decision["choice"],
            "selected_tag": decision["selected_tag"], "fusion_site": decision["fusion_site"],
            "status": "declared_plan_consistency_pass",
            "fusion_junctions": fusion_check,
            "scope": "Does not verify user consent, literature claims or biological neutrality of a tag.",
        })
    missing = required_targets - set(by_target)
    if missing:
        raise ValueError(f"Missing user tag decisions for new expression CDS targets: {sorted(missing)}")
    return checks


def assemble(data, base, allow_test_fixtures=False):
    if data["unresolved_issues"]:
        raise ValueError("Unresolved issues: " + "; ".join(data["unresolved_issues"]))
    source_records, source_meta = {}, {}
    fixture = False
    for source in data["sources"]:
        evidence = source["evidence"]
        if not evidence["reviewed"]:
            raise ValueError(f"{source['id']}: sequence evidence not reviewed.")
        if evidence["evidence_class"] == "test_fixture":
            fixture = True
            if not allow_test_fixtures:
                raise ValueError("Test fixture is not a biological design; explicit test opt-in required.")
        elif evidence["evidence_class"] == "source_record" and not (
            evidence.get("source_url") or evidence.get("accession_version")
        ):
            raise ValueError(f"{source['id']}: source record has neither URL nor accession.")
        path = (Path(base) / source["file"]).resolve()
        rec = read_record(path, source["format"], source.get("record_id"))
        actual_hash = digest(rec.seq)
        if source.get("sequence_sha256") and actual_hash != source["sequence_sha256"]:
            raise ValueError(f"{source['id']}: source sequence hash mismatch.")
        source_records[source["id"]] = rec
        source_meta[source["id"]] = {
            "file": str(path), "file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "sequence_sha256": actual_hash, "evidence": evidence,
        }
    result = SeqRecord(Seq(""), id=data["construct_id"], name=data["construct_id"][:16], description="Computational design, not experimentally validated")
    provenance = []
    for index, block in enumerate(data["assembly"], 1):
        if block["source_id"] not in source_records:
            raise ValueError(f"Unknown source: {block['source_id']}")
        source = source_records[block["source_id"]]
        start, end = block["start"] - 1, block["end"]
        if not 0 <= start < end <= len(source):
            raise ValueError(f"Assembly block {index}: source interval out of bounds.")
        inherit = block.get("inherit_features", True)
        if inherit:
            for f in source.features:
                if f.type == "source":
                    continue
                overlap = any(int(p.end) > start and int(p.start) < end for p in f.location.parts)
                contained = all(start <= int(p.start) and int(p.end) <= end for p in f.location.parts)
                if overlap and not contained:
                    raise ValueError(f"Block {index}: source feature crosses cut; disable inheritance and reannotate explicitly.")
        piece = source[start:end]
        piece.features = [f for f in piece.features if f.type != "source"] if inherit else []
        if block["strand"] == -1:
            piece = piece.reverse_complement(id=True, name=True, description=True)
        for fi, f in enumerate(piece.features, 1):
            f.id = f"block{index}.feature{fi}"
            f.qualifiers["design_id"] = [f.id]
        offset = len(result)
        result = result + piece
        provenance.append({
            "block": index, "source_id": block["source_id"], "source_start": start + 1,
            "source_end": end, "strand": block["strand"], "construct_start": offset + 1,
            "construct_end": offset + len(piece), **source_meta[block["source_id"]],
        })
    if not len(result):
        raise ValueError("Empty assembled construct.")
    result.id = data["construct_id"]
    result.name = data["construct_id"][:16]
    result.description = "TEST FIXTURE ONLY" if fixture else "Computational design; physical identity and function unverified"
    result.annotations = {
        "molecule_type": "DNA", "topology": data["topology"],
        "comment": f"Replication assessment: {data['replication_class']}. {result.description}.",
    }
    result.features.extend(explicit_feature(f, len(result)) for f in data["features"])
    ids = [feature_id(f) for f in result.features]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate final feature IDs.")
    missing = set(data.get("required_feature_ids", [])) - set(ids)
    if missing:
        raise ValueError(f"Missing required features: {sorted(missing)}")
    tag_checks = check_tag_decisions(data, result, source_records)
    translation_checks = check_features(result)
    for region in data.get("protected_regions", []):
        actual = str(location(region, len(result)).extract(result.seq))
        if region.get("sequence") and actual != region["sequence"]:
            raise ValueError(f"Protected region {region['id']} does not match its expected sequence.")
    reference_check = "not_automatically_verified"
    if data["replication_class"] == "exact_plasmid":
        ref_id = data.get("replication_reference_source_id")
        if ref_id not in source_records:
            raise ValueError("exact_plasmid requires a full replication reference source.")
        if not equivalent(result.seq, source_records[ref_id].seq, data["topology"] == "circular"):
            raise ValueError("Assembled sequence is not equivalent to exact-plasmid reference.")
        reference_check = "sequence_equivalence_pass; source identity still requires evidence review"
    qc = {
        "status": "test_fixture_only" if fixture else "sequence_computation_pass",
        "length": len(result), "sequence_sha256": digest(result.seq), "topology": data["topology"],
        "replication_class": data["replication_class"], "reference_check": reference_check,
        "translation_checks": translation_checks, "restriction_sites": sites_report(data, result.seq),
        "tag_decision_checks": tag_checks,
        "experimental_validation": "not_performed", "evo2": "not_performed",
    }
    return result, provenance, qc


def write_tsv(path, columns, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v for k, v in row.items()})


def svg_map(record, path, expression_only=False):
    expression_types = {"CDS", "regulatory", "promoter", "terminator", "intron", "5'UTR", "3'UTR"}
    features = [
        f for f in record.features
        if not expression_only or f.type in expression_types
        or f.qualifiers.get("design_role", [""])[0] in {"target", "tag", "flexible_linker", "signal", "other"}
    ]
    width, height, x0, span = 1200, 160 + 40 * max(1, len(features)), 370, 780
    esc = lambda x: html.escape(str(x), quote=True)
    title = "Expression features" if expression_only else "Full-length annotated map"
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<g font-family="sans-serif" fill="#153448">',
        f'<text x="24" y="28" font-size="18">{esc(record.id)} | {title}</text>',
        f'<text x="24" y="51" font-size="12">{len(record)} bp | {esc(record.annotations.get("topology"))} | COMPUTATIONAL DESIGN, NOT EXPERIMENTALLY VERIFIED</text>',
        f'<line x1="{x0}" y1="82" x2="{x0+span}" y2="82" stroke="#82929f"/>',
        f'<text x="{x0}" y="75" font-size="11">1</text>',
        f'<text x="{x0+span-35}" y="75" font-size="11">{len(record)}</text>',
    ]
    if not features:
        lines.append('<text x="24" y="115" font-size="13">No expression features annotated; inspect the full map and evidence.</text>')
    for i, f in enumerate(features):
        y = 110 + 40 * i
        label = f.qualifiers.get("label", f.qualifiers.get("gene", [feature_id(f)]))[0]
        text = f"{i+1}. {label} [{f.type}]"
        short = text if len(text) < 47 else text[:44] + "..."
        coords = ",".join(f"{int(p.start)+1}..{int(p.end)}" for p in f.location.parts)
        color = "#167d9a" if f.type == "CDS" else "#9b6725" if f.type == "intron" else "#53688b"
        lines += [
            f'<text x="24" y="{y+4}" font-size="12"><title>{esc(text)}</title>{esc(short)}</text>',
            f'<text x="24" y="{y+18}" font-size="10">{esc(coords)} | strand {f.location.strand}</text>',
        ]
        for pi, p in enumerate(f.location.parts):
            left = x0 + int(p.start) / len(record) * span
            right = x0 + int(p.end) / len(record) * span
            extent = max(1, right - left)
            lines.append(f'<rect x="{left:.3f}" y="{y-7}" width="{extent:.3f}" height="14" fill="{color}"><title>{esc(text)} | {esc(coords)}</title></rect>')
            if pi == len(f.location.parts) - 1 and p.strand in {1, -1}:
                tip = left if p.strand == -1 else right
                tail = tip + min(6, extent) if p.strand == -1 else tip - min(6, extent)
                lines.append(f'<polygon points="{tip:.3f},{y} {tail:.3f},{y-7} {tail:.3f},{y+7}" fill="#163448"/>')
    lines.append("</g></svg>")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def build(manifest_path, out, allow_test_fixtures=False):
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be new or empty.")
    data = manifest(manifest_path)
    record, provenance, qc = assemble(data, Path(manifest_path).resolve().parent, allow_test_fixtures)
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".vector-export-", dir=out.parent) as temporary:
        stage = Path(temporary)
        SeqIO.write(record, stage / "construct.gb", "genbank")
        SeqIO.write(record, stage / "construct.fasta", "fasta")
        reread = read_record(stage / "construct.gb", "genbank")
        fa = read_record(stage / "construct.fasta", "fasta")
        if reread.seq != record.seq or fa.seq != record.seq or reread.annotations["topology"] != data["topology"]:
            raise ValueError("Sequence/topology round-trip failed.")
        if len(reread.features) != len(record.features):
            raise ValueError("Feature count round-trip failed.")
        for original, roundtrip in zip(record.features, reread.features):
            if original.location != roundtrip.location or original.extract(record.seq) != roundtrip.extract(reread.seq):
                raise ValueError("Feature coordinate/extraction round-trip failed.")
        check_features(reread)
        qc["file_roundtrip"] = "pass"
        rows = [{
            "feature_id": feature_id(f), "type": f.type,
            "segments_1based_inclusive": [[int(p.start)+1, int(p.end)] for p in f.location.parts],
            "strand": f.location.strand, "qualifiers": f.qualifiers,
        } for f in reread.features]
        write_tsv(stage / "features.tsv", ["feature_id", "type", "segments_1based_inclusive", "strand", "qualifiers"], rows)
        write_tsv(stage / "sources.tsv", ["block", "source_id", "file", "source_start", "source_end", "strand", "construct_start", "construct_end", "file_sha256", "sequence_sha256", "evidence"], provenance)
        svg_map(reread, stage / "plasmid_map.svg")
        svg_map(reread, stage / "expression_map.svg", True)
        dump_json(stage / "qc_report.json", qc)
        (stage / "README.md").write_text(
            f"# {record.id}\n\n{record.description}. All coordinates in features.tsv are 1-based inclusive; joined segments are in biological order.\n"
            "These files describe a computed sequence, not a sequenced physical sample. See qc_report.json and source evidence before use.\n",
            encoding="utf-8",
        )
        if out.exists():
            out.rmdir()
        out.mkdir()
        for item in stage.iterdir():
            item.rename(out / item.name)
    return qc


def compare_sequences(data, baseline, candidate, constraint_features=None):
    baseline, candidate = dna(baseline), dna(candidate)
    report = {
        "baseline_sha256": digest(baseline), "candidate_sha256": digest(candidate),
        "baseline_length": len(baseline), "candidate_length": len(candidate),
        "status": "fail", "issues": [], "checks": [],
    }
    if len(baseline) != len(candidate):
        report["issues"].append("Different lengths: new coordinate mapping and candidate annotation required.")
        report["comparability"] = "incomparable_without_reannotation"
        return report
    report["comparability"] = "equal_length_coordinate_checks_only"
    report["differences"] = [
        {"pos": i + 1, "ref": a, "alt": b}
        for i, (a, b) in enumerate(zip(baseline, candidate)) if a != b
    ]
    for region in data.get("protected_regions", []):
        loc = location(region, len(baseline))
        ref, alt = str(loc.extract(Seq(baseline))), str(loc.extract(Seq(candidate)))
        if ref != alt:
            report["issues"].append(f"Protected region changed: {region['id']}")
        if region.get("sequence") and ref != region["sequence"]:
            report["issues"].append(f"Baseline protected-region expectation mismatch: {region['id']}")
    features = constraint_features if constraint_features is not None else [
        explicit_feature(spec, len(baseline)) for spec in data["features"]
    ]
    for f in features:
        if f.type != "CDS":
            continue
        try:
            expected_values = f.qualifiers.get("translation")
            if not expected_values:
                raise ValueError(f"{feature_id(f)}: no expected protein for candidate constraints.")
            expected = expected_values[0].replace(" ", "").replace("\n", "").rstrip("*")
            if protein(f, baseline) != expected or protein(f, candidate) != expected:
                report["issues"].append(f"Protein identity changed/mismatched: {feature_id(f)}")
        except Exception as exc:
            report["issues"].append(str(exc))
    try:
        report["restriction_sites"] = sites_report(data, candidate)
    except ValueError as exc:
        report["issues"].append(str(exc))
    report["baseline_gc"] = (baseline.count("G") + baseline.count("C")) / len(baseline)
    report["candidate_gc"] = (candidate.count("G") + candidate.count("C")) / len(candidate)
    report["status"] = "pass" if not report["issues"] else "fail"
    report["scope"] = "Declared constraints and provided feature constraints only; no unannotated functions are inferred."
    return report


def apply_edits(data, baseline, edits, constraint_features=None):
    baseline = dna(baseline)
    touched = set()
    normalized = []
    for edit in edits:
        if set(edit) != {"pos", "ref", "alt"} or type(edit["pos"]) is not int:
            raise ValueError("Each edit requires exactly integer pos, ref and alt.")
        start, ref, alt = edit["pos"] - 1, dna(edit["ref"]), dna(edit["alt"])
        if len(ref) != len(alt):
            raise ValueError("Only equal-length replacements are supported.")
        if start < 0 or start + len(ref) > len(baseline) or baseline[start:start+len(ref)] != ref:
            raise ValueError("Edit REF mismatch or out of bounds.")
        positions = set(range(start, start + len(ref)))
        if touched & positions:
            raise ValueError("Overlapping edits are not allowed.")
        touched.update(positions)
        normalized.append((start, ref, alt))
    candidate = baseline
    for start, ref, alt in sorted(normalized, reverse=True):
        candidate = candidate[:start] + alt + candidate[start+len(ref):]
    report = compare_sequences(data, baseline, candidate, constraint_features)
    if report["status"] != "pass":
        raise ValueError("; ".join(report["issues"]))
    return candidate, report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("build")
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--allow-test-fixtures", action="store_true")
    for name in ("compare", "apply-edits"):
        p = commands.add_parser(name)
        p.add_argument("--manifest", type=Path, required=True)
        p.add_argument("--baseline", type=Path, required=True)
        p.add_argument("--out", type=Path, required=True)
        p.add_argument("--candidate" if name == "compare" else "--edits", type=Path, required=True)
        p.add_argument("--allow-test-fixtures", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "build":
            qc = build(args.manifest, args.out, args.allow_test_fixtures)
            print(json.dumps(qc, ensure_ascii=False))
        else:
            if args.out.exists() or (
                args.command == "apply-edits" and args.out.with_suffix(args.out.suffix + ".qc.json").exists()
            ):
                raise ValueError("Output already exists; use a new versioned filename.")
            data = manifest(args.manifest)
            baseline = read_record(args.baseline)
            original, _, _ = assemble(data, args.manifest.resolve().parent, args.allow_test_fixtures)
            if baseline.seq != original.seq:
                raise ValueError("Baseline does not match the frozen manifest assembly.")
            if args.command == "compare":
                candidate = read_record(args.candidate)
                report = compare_sequences(data, baseline.seq, candidate.seq, original.features)
                args.out.parent.mkdir(parents=True, exist_ok=True)
                dump_json(args.out, report)
                return 0 if report["status"] == "pass" else 1
            edits = json.loads(args.edits.read_text(encoding="utf-8"))
            candidate, report = apply_edits(data, baseline.seq, edits, original.features)
            args.out.parent.mkdir(parents=True, exist_ok=True)
            SeqIO.write(SeqRecord(Seq(candidate), id=baseline.id + "_candidate", description="QC-filtered equal-length candidate; not experimentally validated"), args.out, "fasta")
            dump_json(args.out.with_suffix(args.out.suffix + ".qc.json"), report)
        return 0
    except Exception as exc:
        if args.command == "build" and (not args.out.exists() or not any(args.out.iterdir())):
            args.out.mkdir(parents=True, exist_ok=True)
            dump_json(args.out / "qc_report.json", {"status": "fail", "issues": [str(exc)], "final_sequence_exported": False})
        print(json.dumps({"status": "fail", "error": str(exc)}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
