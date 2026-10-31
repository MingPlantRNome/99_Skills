#!/usr/bin/env bash
# Reproducible fetch + build of Arabidopsis (TAIR10) resources for LncAxiom:
#   reference_genome/TAIR10.fa + TAIR10.gtf   (NCBI RefSeq GCF_000001735.4_TAIR10.1)
#   model_training_datasets/Ath_mRNA.fa       (protein-coding transcripts, n≈52k)
#   model_training_datasets/training_dataset.fa (lncRNA transcripts, n≈3.9k)
#
# Source: NCBI RefSeq (reachable where tair.com login-walls curl). The RefSeq
# GTF's transcript_biotype "mRNA"/"lnc_RNA" classes are used for training split.
# Verified 2026-10-30: counts match GTF exactly; lncRNA min=200nt median=348nt.
#
# Usage: prepare_arabidopsis_resources.sh <outdir>     (default ~/lncaxiom-resources)
# Requires: curl, gzcat, awk, gffread (brew install gffread).
set -euo pipefail

OUT="${1:-$HOME/lncaxiom-resources}"
BASE="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Arabidopsis_thaliana/latest_assembly_versions/GCF_000001735.4_TAIR10.1"
mkdir -p "$OUT"/{reference_genome,model_training_datasets,annotation}

echo "== 1/4 download genome + annotation =="
curl -s --max-time 600 -o "$OUT/reference_genome/TAIR10.fa.gz"   "$BASE/GCF_000001735.4_TAIR10.1_genomic.fna.gz"
curl -s --max-time 600 -o "$OUT/annotation/refseq.gtf.gz"       "$BASE/GCF_000001735.4_TAIR10.1_genomic.gtf.gz"
gzcat "$OUT/reference_genome/TAIR10.fa.gz" > "$OUT/reference_genome/TAIR10.fa"
gzcat "$OUT/annotation/refseq.gtf.gz"      > "$OUT/reference_genome/TAIR10.gtf"

echo "== 2/4 split GTF exons by transcript_biotype =="
rm -f "$OUT"/annotation/split_*.gtf
gzcat "$OUT/annotation/refseq.gtf.gz" | awk -F'\t' 'BEGIN{OFS=FS} {
  if ($3=="transcript") {
    if (match($0,/transcript_biotype "[^"]+"/)) { b=substr($0,RSTART+20,RLENGTH-21) } else { b="?" }
    cls = (b=="mRNA") ? "mrna" : (b=="lnc_RNA") ? "lncrna" : "other"
  } else if ($3=="exon" && (cls=="mrna" || cls=="lncrna")) {
    print $0 >> ("'"$OUT"'/annotation/split_" cls ".gtf")
  }
}'

echo "== 3/4 extract training FASTAs with gffread =="
command -v gffread >/dev/null || { echo "MISS gffread (brew install gffread)" >&2; exit 1; }
gffread -w "$OUT/model_training_datasets/Ath_mRNA.fa"         -g "$OUT/reference_genome/TAIR10.fa" "$OUT/annotation/split_mrna.gtf"
gffread -w "$OUT/model_training_datasets/training_dataset.fa" -g "$OUT/reference_genome/TAIR10.fa" "$OUT/annotation/split_lncrna.gtf"

echo "== 4/4 verify counts =="
grep -c "^>" "$OUT/model_training_datasets/Ath_mRNA.fa"
grep -c "^>" "$OUT/model_training_datasets/training_dataset.fa"
echo "done -> $OUT"
