---
name: lncaxiom-lncrna-pipeline
description: Run the LncAxiom Snakemake workflow for genome-wide lncRNA identification, annotation, and functional prioritization from RNA-seq (HISAT2/StringTie/FEELnc/CPC2/DESeq2/psRobot). Use when a transcriptome-based lncRNA discovery pipeline with modular stages is needed. Requires snakemake >=7.8 and conda/mamba; Linux or macOS.
tool_type: workflow
primary_tool: LncAxiom (snakemake)
---

# LncAxiom lncRNA Pipeline (WRAP)

Thin ResearchOS wrapper around the upstream **LncAxiom** Snakemake workflow. This
skill contains orchestration knowledge only — the workflow code lives upstream and
is fetched at run time. Do not copy upstream rules into projects.

> **Local patches (2026-10-31)** — `scripts/install_resources.sh` and
> `scripts/lncaxiom_run.sh` in this bundle carry two fixes, verified on a fresh
> Apple Silicon Mac (mlys-M5). Details in `PATCHES.md`:
> 1. `install_resources.sh`: PLncDB FASTA check now tolerates the Zenodo layout
>    (`homology_databases/PLncDB/A.thaliana/…`) in addition to the flat one.
> 2. `lncaxiom_run.sh`: snakemake is located by probing real `/Users/<user>`
>    (CSSwitch sandbox virtualises `$HOME`) and common conda roots, not only the
>    bare `command -v` / `$HOME/miniforge3`. The resolved binary runs the
>    entire dry-run/execute path. Require `conda` on PATH for `--use-conda`.
> Note: upstream BLAST-novelty module has a DAG gap (rule `blast_against_db`
> requires `results/feelnc/candidates.fa`, which no rule produces) — keep
> `run_blast_novelty: false` unless the upstream is fixed.

## Upstream provenance

| Field | Value |
|---|---|
| Repository | `numanmian/LncAxiom-` (https://github.com/numanmian/LncAxiom-) |
| License | MIT (LICENSE verified in-tree) |
| Requires | snakemake >=7.8, conda/mamba, git; Linux/macOS |
| Imported from | user attachment `LncAxiom--main.zip` (= repo `main` branch, 2026-10) |
| Registry | `registry/upstreams.jsonl` → `github:numanmian/LncAxiom-`, ADOPT+WRAP |

## When to use

- Genome-wide **transcriptome-based** lncRNA discovery from RNA-seq: preprocessing →
  HISAT2 alignment → StringTie assembly → merge/quantify → candidate filtering →
  coding-potential prediction (FEELnc + CPC2) → classification → DESeq2 differential
  expression → (optional) psRobot eTM / BLAST novelty / functional enrichment.
- You need a **modular** run: module switches in `config.yaml` let you re-run only
  downstream stages without repeating upstream steps.

## When NOT to use

- Genome-assembly ncRNA scanning (tRNA/rRNA/snoRNA via Rfam/Infernal) — use
  `bio-genome-annotation-ncrna-annotation` instead (different task, reviewed skill).
- Standalone differential expression of known genes — use `bio-de-deseq2-basics`.
- You only have a genome and no RNA-seq: this pipeline is transcriptome-driven.

## Input contract

All configured in upstream `config/config.yaml` + `config/samples.tsv`:

| Input | Required | Notes |
|---|---|---|
| `samples.tsv` | yes | Columns: `sample, layout, condition, sra`; layout SINGLE/PAIR; condition levels must match `contrast` in config |
| Genome FASTA + GTF | yes | Default points at TAIR10 (Arabidopsis) |
| RNA-seq reads | yes | `input_method: sra` (auto-download) or `local` (`{sample}.fastq.gz` in `local_fastq_dir`) |
| FEELnc training FASTAs | yes for prediction | **NOT shipped upstream** — zip `resources/` only contains miRNAs.fa; you must supply both files. Default paths expect Arabidopsis (`Ath_mRNA.fa` + `training_dataset.fa`) — swap for other species, coding-potential results are species-dependent. For Arabidopsis run `scripts/prepare_arabidopsis_resources.sh <outdir>` — downloads TAIR10 (NCBI RefSeq GCF_000001735.4) and builds both training FASTAs from the `transcript_biotype` mRNA/lnc_RNA classes (verified: 52,179 mRNA / 3,878 lncRNA, lncRNA min 200 nt). |
| miRNA FASTA | only if psRobot eTM on | `run_psrobot_etm: true` |
| BLAST dbs (CANTATA/NONCODE/PLncDB/GreeNc) | only if novelty on | see "Homology databases" below — CANTATAdb 3.0 + PLncDB acquired and BLAST-indexed locally; **GreeNC unavailable** (down since ≥2024-03, confirmed by the CANTATAdb 3.0 paper); NONCODE not yet fetched. Edit `blast_databases` in config to list only existing DBs. |

## Distribution bundle (portable install)

**The skill bundle ships WITHOUT large data files** (installers such as CSSwitch
reject bundles containing files > 4 MiB). Genome, training FASTAs and BLAST
indexes live in a separate **resources directory** that you pass as a parameter:

```bash
# 1. obtain a resources directory (pre-built, or regenerate with
#    scripts/prepare_arabidopsis_resources.sh), e.g.:
RES=~/lncaxiom-resources            # classic layout
RES=~/.CSSwitch/sandbox/home/lncaxiom-resources   # inside the CSSwitch sandbox

# 2. verify it and get the config block:
bash scripts/install_resources.sh "$RES"

# 3. run — the adapter wires config paths from --resources automatically:
bash scripts/lncaxiom_run.sh <workdir> --resources "$RES" \
    [--from-bam /path/to/bam_dir] [--cores N]
# (paste the printed blast_databases block into config.yaml on first run)
```

The git repository copy of this skill intentionally holds no data at all —
`resources/` only ever exists in standalone bundles or sandbox deployments.

## Homology databases (novelty BLAST) status — 2026-10-30

Pre-built BLAST DBs (makeblastdb format) live in the resources directory
(`resources/homology_databases/` of a bundle/sandbox deployment).
Built 2026-10-30 with NCBI BLAST+ 2.17.0 standalone (fetched from
`ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/`; rebuilding only
requires any blast+ >= 2.10):

| DB | Content | Count | Source |
|---|---|---|---|
| `CANTATAdb/CANTATAdb.*` | A.thaliana lncRNAs, CANTATAdb 3.0 | 6,775 | `rhesus.amu.edu.pl/CANTATA/DOWNLOADS/Arabidopsis_thaliana_lncRNAs.fasta` (count matches Szcześniak & Wanowska 2024, Plant Cell Physiol 65:1486, DOI 10.1093/pcp/pcae081) |
| `PLncDB/PLncDB.*` | A.thaliana lncRNAs, PLncDB core set | 1,682 | Zenodo backup of PLncDB (`zenodo.org/record/4017591`), file `A.thaliana.tar.gz` → `A.tha.lncRNA-PLncDB.fa`. Optional larger sets in same tarball: Liu 7,887 / RNAcentral 3,886 / EVLncRNAs 144 |
| ~~GreeNc~~ | — | — | **unavailable**: host down (500/timeout), old host DNS dead; CANTATAdb 3.0 paper states "not accessible at the time of writing (March 2024)". No mirror found. |
| NONCODE | — | — | not fetched yet (add when enabling novelty for non-Arabidopsis work) |

To enable the novelty module: point `blast_databases` at these two DBs only, delete the GreeNc/NONCODE entries from config, set `run_blast_novelty: true`.

## Safe execution workflow

```bash
# Thin adapter: fetch upstream, preflight-check inputs, dry-run, then execute.
bash scripts/lncaxiom_run.sh <workdir> [--cores N] [--skip-dry-run]
```

Or manually inside a fresh clone of `numanmian/LncAxiom-`:

```bash
snakemake -n --use-conda          # mandatory dry run first
snakemake --use-conda --cores 8
```

Rules of engagement:

- **Never auto-install** missing genomes/DBs/training data — stop and report
  (`scripts/env_capture.py check --tool snakemake` confirms the runner only).
- Per-rule conda envs come from upstream `envs/*.yaml` via `--use-conda`; first
  run builds them (slow, several GB).
- Raw FASTQ / SRA data are **P2** — this pipeline runs fully locally; do not paste
  paths or contents into cloud chats.
- Record the run: `run_ledger.py init` + `env_capture.py capture --run RUN-...`.
- HPC/remote submission requires explicit user approval (POLICY.md).

## Starting from external BAMs (skip download/QC/mapping)

If mapping was done by collaborators, start directly at **assembly** — no upstream
code changes needed, because module includes are config-gated in the Snakefile
and every downstream rule consumes only
`results/aligned_bam/{sample}.sorted.bam`(+`.bai`):

```bash
bash scripts/lncaxiom_run.sh <workdir> --from-bam /path/to/bam_dir [--cores N]
```

Prerequisites for `/path/to/bam_dir`:

- One `{sample}.sorted.bam` **and** `{sample}.sorted.bam.bai` per sample, where
  `{sample}` matches the `sample` column of `config/samples.tsv` exactly.
- BAMs must be **coordinate-sorted and indexed**; the adapter runs
  `samtools quickcheck` on all of them and hard-fails on corruption.
- **Reference genome must match the pipeline config**: contig names in the BAM
  header must equal those in the configured `genome_gtf`/`genome_fasta`
  (TAIR10 by default). Mismatched naming (e.g. `Chr1` vs `chr1` vs `1`) is the
  #1 cause of empty downstream GTFs — verify with
  `samtools view -H x.bam | grep '^@SQ'` before running.

What `--from-bam` does, in order:

1. Sets `run_download: false` and `run_preprocessing: false` in the clone's
   `config/config.yaml` (mapping + fastqc/fastp rules are excluded from the DAG;
   the `rule all` targets for those modules disappear automatically).
2. Verifies every sample's `.sorted.bam`/`.bai` exists (fails with a list if not).
3. Symlinks them into `results/aligned_bam/` (never copies or modifies).
4. `samtools quickcheck` all BAMs (warn-and-skip only if samtools is absent).
5. Dry-runs (or executes with `--skip-dry-run`) — snakemake builds the DAG
   starting at StringTie assembly.

Manual equivalent: flip the two switches yourself and place BAMs at
`results/aligned_bam/{sample}.sorted.bam[.bai]`, then `snakemake -n --use-conda`.

Caveats when entering at BAM:

- No fastqc/fastp/multiqc evidence exists for your samples — obtain QC summary
  from whoever did the mapping, or run `samtools flagstat`/`idxstats` yourself;
  record alignment rates as diagnostic figures before biological conclusions.
- `input_method` in config becomes irrelevant (it only gates the download rule).
- Downstream stages (quantification → prediction → FEELnc → DESeq2) run
  unchanged; re-runs of later modules alone use the module switches below.

## Module switches (config.yaml)

`run_download/preprocessing/assembly/quantification/lncrna_prediction`,
`run_feelnc_filter/codpot/classifier`, `run_psrobot_etm`, `run_blast_novelty`,
`run_enrichment`. Defaults run the full Arabidopsis demo path with novelty/psRobot
off. DESeq2 contrast (`factor/case/control`) **must be edited** to your design.

## Expected outputs (QC before biology)

- `results/merged.gtf`, `results/gene_counts.csv` / `transcript_counts.csv`
- `results/lncrna_candidates.gtf`, FEELnc high-confidence GTF
- DESeq2 results + volcano (padj<0.05, |LFC|>1 defaults)

Diagnostic figures before concluding: per-sample alignment/assignment rates,
count-depth distributions, candidate-count sanity vs. literature scale for the
species, volcano/MA plot. Bind lncRNA-count claims to a Figure record
(`figure_tool.py`) before writing them into a brief.

## Common failure modes (from upstream + snakemake practice)

1. FEELnc codpot fails with missing `train_mrna_fasta`/`train_lncrna_fasta` —
   the training FASTAs are **not shipped** with the workflow; prepare them per
   the input contract before enabling prediction.
2. `Can't locate Bio/DB/SeqFeature/Store.pm` — FEELnc classifier missing
   BioPerl module; install inside the active conda env, re-run the failed rule.
3. `No executable svm-scale on CPC2 path!` — CPC2/libsvm not compiled in env;
   reinstall CPC2 in its rule env.
4. Empty candidate GTF — check StringTie merge succeeded and GTF chromosome
   names match the genome FASTA (`grep '^>'` vs GTF seqnames).
5. Contrast mismatch — DESeq2 fails when `case`/`control` levels are absent
   from `samples.tsv` `condition`.
6. First-run conda build failures — pin or re-run; `--conda-frontend mamba`
   speeds env creation.

## Validation

- `snakemake -n` completes with the expected target count (dry-run DAG sane).
- Upstream ships `test_data/` for a minimal smoke path.
- After a real run: `env_capture.py capture` + record upstream commit in the Run
  ledger (`upstreams` field) so results are reproducible.
