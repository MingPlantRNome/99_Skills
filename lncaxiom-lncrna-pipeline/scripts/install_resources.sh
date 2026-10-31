#!/usr/bin/env bash
# Verify a LncAxiom resources directory and print the config.yaml block for it.
# The skill bundle ships WITHOUT large data — data lives in a separate resources
# directory whose path is the parameter.
#
# Usage:
#   install_resources.sh [resources_dir]
#     default ./resources next to this script (classic full-bundle layout)
#   CSSwitch sandbox:  install_resources.sh ~/.CSSwitch/sandbox/home/lncaxiom-resources
set -euo pipefail
RES="${1:-$(cd "$(dirname "$0")/.." && pwd)/resources}"
[ -d "$RES" ] || { echo "resources dir not found: $RES" >&2; exit 1; }

echo "== 1/3 decompress large files (if still gzipped) =="
for f in "$RES/reference_genome/TAIR10.fa" "$RES/reference_genome/TAIR10.gtf" \
         "$RES/model_training_datasets/Ath_mRNA.fa"; do
  [ -f "$f" ] || gunzip -k "$f.gz"
done

echo "== 2/3 verify counts =="
check() { n=$(grep -c '^>' "$1"); echo "  $(basename "$1"): $n seqs (expected $2)"; [ "$n" -eq "$2" ] || { echo "COUNT MISMATCH: $1" >&2; exit 1; }; }
check "$RES/model_training_datasets/Ath_mRNA.fa" 52179
check "$RES/model_training_datasets/training_dataset.fa" 3878
check "$RES/homology_databases/CANTATAdb/Arabidopsis_thaliana_lncRNAs.fasta" 6775
# PLncDB fasta lives either flat in PLncDB/ (classic layout) or under the
# Zenodo-extracted A.thaliana/ subdir. Probe both before checking.
PLNCDB_FA="$RES/homology_databases/PLncDB/A.tha.lncRNA-PLncDB.fa"
[ -f "$PLNCDB_FA" ] || PLNCDB_FA="$RES/homology_databases/PLncDB/A.thaliana/A.tha.lncRNA-PLncDB.fa"
check "$PLNCDB_FA" 1682
for db in CANTATAdb/CANTATAdb PLncDB/PLncDB; do
  [ -f "$RES/homology_databases/$db.nin" ] || { echo "MISS BLAST index: $db" >&2; exit 1; }
done
echo "  BLAST indexes present (CANTATAdb, PLncDB)"

echo "== 3/3 config.yaml block =="
cat <<EOF
# paste into <workdir>/config/config.yaml:
genome_fasta:         "$RES/reference_genome/TAIR10.fa"
genome_gtf:           "$RES/reference_genome/TAIR10.gtf"
train_lncrna_fasta:   "$RES/model_training_datasets/training_dataset.fa"
train_mrna_fasta:     "$RES/model_training_datasets/Ath_mRNA.fa"
psrobot_mirna_fasta:  "$RES/microrna_sequences/miRNAs.fa"
run_blast_novelty:    true
blast_databases:
  cantatadb: "$RES/homology_databases/CANTATAdb/CANTATAdb"
  plncdb:    "$RES/homology_databases/PLncDB/PLncDB"
# NOTE: GreeNC entry intentionally absent — database offline since 2024-03.
EOF
echo "done. resources at: $RES"
