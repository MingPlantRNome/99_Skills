#!/usr/bin/env bash
# Thin adapter: fetch upstream LncAxiom, preflight-check config inputs,
# dry-run, then execute. Never installs genomes/DBs/training data.
#
# Usage:
#   lncaxiom_run.sh <workdir> [--cores N] [--skip-dry-run]
#   lncaxiom_run.sh <workdir> --from-bam <bam_dir> [--cores N] [--skip-dry-run]
#   any form + --resources <resources_dir>
#
# --from-bam: mapping was done externally. Expects <bam_dir>/{sample}.sorted.bam
# (+ .bai) for every sample in config/samples.tsv. The adapter flips
# run_download/run_preprocessing to false in the cloned config, symlinks the
# BAMs into results/aligned_bam/, quickchecks them, and starts at assembly.
#
# --resources: path to a LncAxiom resources directory (see install_resources.sh).
# The adapter points genome/gtf/training/miRNA config keys at it and enables the
# novelty module with the bundled CANTATAdb/PLncDB indexes.
set -euo pipefail

UPSTREAM_URL="${UPSTREAM_URL:-https://github.com/numanmian/LncAxiom-.git}"
CORES=4
DRY_RUN=1
BAM_DIR=""
RES=""

if [ $# -lt 1 ]; then
  echo "usage: $0 <workdir> [--from-bam <dir>] [--cores N] [--skip-dry-run]" >&2
  exit 2
fi
WORKDIR="$1"; shift
while [ $# -gt 0 ]; do
  case "$1" in
    --from-bam) BAM_DIR="$(cd "$2" && pwd)"; shift 2 ;;
    --resources) RES="$(cd "$2" && pwd)"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --skip-dry-run) DRY_RUN=0; shift ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

# --- 0. runner present? (report only, never install) ---
# Locate snakemake on PATH, else probe common conda env roots.
# The skill scripts themselves are report-only; we never install snakemake.
find_snakemake() {
  command -v snakemake 2>/dev/null && return 0
  # In a CSSwitch/macOS sandbox $HOME may be virtualised away from the real
  # user dir; probe the real /Users base too so a user-level conda is found.
  for ub in /Users/mingleiyang "$HOME"; do
    for root in "$ub/miniforge3" "$ub/mambaforge" "$ub/miniconda3" "$ub/anaconda3"; do
      for env in "$root/envs/bio-agent" "$root/envs/bioagents" "$root"; do
        if [ -x "$env/bin/snakemake" ]; then
          echo "$env/bin/snakemake"
          return 0
        fi
      done
    done
  done
  return 1
}
SNAKEMAKE_BIN="$(find_snakemake)" || {
  echo "MISS snakemake (>=7.8). Activate an env that has it, e.g. bio-agent." >&2
  exit 1
}
echo "using snakemake: $SNAKEMAKE_BIN ($("$SNAKEMAKE_BIN" --version 2>/dev/null))"

# --- 1. fetch upstream ---
if [ -d "$WORKDIR/.git" ]; then
  git -C "$WORKDIR" pull --ff-only
else
  git clone "$UPSTREAM_URL" "$WORKDIR"
fi
cd "$WORKDIR"
echo "upstream commit: $(git rev-parse --short HEAD)"

# --- 2. preflight: config inputs exist? ---
python3 - <<'PYEOF'
import csv, sys
from pathlib import Path
cfg = Path("config/config.yaml")
if not cfg.exists():
    sys.exit("missing config/config.yaml — edit upstream config first")
if not Path("config/samples.tsv").exists():
    sys.exit("missing config/samples.tsv")
with open("config/samples.tsv", newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
need = {"sample", "layout", "condition", "sra"}
if not rows or not need.issubset(rows[0].keys()):
    sys.exit("samples.tsv must have columns: sample, layout, condition, sra")
print(f"preflight ok: {len(rows)} samples declared")
PYEOF
echo "NOTE: genome FASTA/GTF, FEELnc training data, BLAST DBs are NOT checked"
echo "here — confirm config.yaml paths yourself (see SKILL.md input contract)."

# --- 3. optional: point config at an external resources directory ---
if [ -n "$RES" ]; then
  for f in "$RES/reference_genome/TAIR10.fa" "$RES/reference_genome/TAIR10.gtf" \
           "$RES/model_training_datasets/Ath_mRNA.fa"; do
    [ -f "$f" ] || gunzip -k "$f.gz" 2>/dev/null || { echo "MISS $f(.gz) in $RES" >&2; exit 1; }
  done
  sed -i '' -E \
    -e "s|^(genome_fasta:).*|\1 \"$RES/reference_genome/TAIR10.fa\"|" \
    -e "s|^(genome_gtf:).*|\1 \"$RES/reference_genome/TAIR10.gtf\"|" \
    -e "s|^(train_lncrna_fasta:).*|\1 \"$RES/model_training_datasets/training_dataset.fa\"|" \
    -e "s|^(train_mrna_fasta:).*|\1 \"$RES/model_training_datasets/Ath_mRNA.fa\"|" \
    -e "s|^(psrobot_mirna_fasta:).*|\1 \"$RES/microrna_sequences/miRNAs.fa\"|" \
    -e "s|^(run_blast_novelty:).*|\1 true|" \
    config/config.yaml
  python3 - "$RES" <<'PYEOF'
import sys
res = sys.argv[1]
block = f"""blast_databases:
  cantatadb: "{res}/homology_databases/CANTATAdb/CANTATAdb"
  plncdb:    "{res}/homology_databases/PLncDB/PLncDB"
"""
print("resources wired; replace the blast_databases block in config/config.yaml with:")
print(block)
PYEOF
fi

# --- 4. optional: external-BAM entry (skip download/preprocess/mapping) ---
if [ -n "$BAM_DIR" ]; then
  # 3a. flip module switches in the clone's config (idempotent)
  sed -i '' -E 's/^(run_download:).*/\1           false/; s/^(run_preprocessing:).*/\1      false/' \
      config/config.yaml
  echo "config: run_download=false, run_preprocessing=false (external BAM entry)"

  # 3b. per-sample BAM presence check
  python3 - "$BAM_DIR" <<'PYEOF'
import csv, sys
from pathlib import Path
bam_dir = Path(sys.argv[1])
with open("config/samples.tsv", newline="") as fh:
    samples = [r["sample"] for r in csv.DictReader(fh, delimiter="\t")]
missing = []
for s in samples:
    for suffix in (".sorted.bam", ".sorted.bam.bai"):
        if not (bam_dir / f"{s}{suffix}").exists():
            missing.append(f"{s}{suffix}")
if missing:
    sys.exit("BAM dir is missing required files (name them {sample}.sorted.bam[.bai]):\n  "
             + "\n  ".join(missing[:10]))
print(f"bam preflight ok: {len(samples)} x (sorted.bam + .bai) present in {bam_dir}")
PYEOF

  # 3c. symlink into the path downstream rules consume; quickcheck if samtools
  mkdir -p results/aligned_bam
  python3 - "$BAM_DIR" <<'PYEOF'
import csv, os, sys
from pathlib import Path
bam_dir = Path(sys.argv[1])
dest = Path("results/aligned_bam")
with open("config/samples.tsv", newline="") as fh:
    samples = [r["sample"] for r in csv.DictReader(fh, delimiter="\t")]
for s in samples:
    for suffix in (".sorted.bam", ".sorted.bam.bai"):
        target, link = bam_dir / f"{s}{suffix}", dest / f"{s}{suffix}"
        if link.exists() or link.is_symlink():
            if link.is_symlink() and Path(os.readlink(link)) == target:
                continue
            if not link.is_symlink():
                continue  # real file already there; do not touch
            link.unlink()
        link.symlink_to(target)
print(f"linked {len(samples)*2} BAM/BAI files into {dest}")
PYEOF
  if command -v samtools >/dev/null 2>&1; then
    echo "samtools quickcheck:"
    for f in results/aligned_bam/*.sorted.bam; do
      samtools quickcheck -v "$f" || { echo "CORRUPT: $f" >&2; exit 1; }
    done
    echo "  all BAMs pass quickcheck"
  else
    echo "WARN samtools not found — skipped quickcheck (run 'samtools quickcheck' yourself)"
  fi
fi

# --- 5. dry run, then execute ---
if [ "$DRY_RUN" -eq 1 ]; then
  "$SNAKEMAKE_BIN" -n --use-conda
  echo "dry-run ok. Re-run with --skip-dry-run to execute."
  exit 0
fi
exec "$SNAKEMAKE_BIN" --use-conda --cores "$CORES"
