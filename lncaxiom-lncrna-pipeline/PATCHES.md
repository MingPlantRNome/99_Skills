# LncAxiom skill — local patches

Two fixes to `scripts/install_resources.sh` and `scripts/lncaxiom_run.sh`,
verified 2026-10-31 on a fresh Apple Silicon Mac (macOS, mlys-M5). They were
discovered while running the "new-machine 3-step" flow against
`/Users/mingleiyang/lncaxiom-resources` under a CSSwitch sandbox.

## Fix 1 — install_resources.sh: PLncDB FASTA layout
- **Symptom**: the 4-count verification fails at the PLncDB entry with
  `grep: <...>/PLncDB/A.tha.lncRNA-PLncDB.fa: No such file or directory`
  (exit 2) even though the file exists and has the expected 1682 sequences.
- **Cause**: the script hard-coded `homology_databases/PLncDB/A.tha.lncRNA-PLncDB.fa`,
  but Zenodo-backed resource bundles put it under
  `homology_databases/PLncDB/A.thaliana/A.tha.lncRNA-PLncDB.fa` (one extra
  `A.thaliana/` level).
- **Fix**: probe the flat path first, fall back to the `A.thaliana/` subdir
  before calling `check`.

## Fix 2 — lncaxiom_run.sh: snakemake discovery
- **Symptom**: the adapter exits `MISS snakemake (>=7.8)` right away, even on
  machines that have snakemake in a conda env.
- **Cause**: it relied on `command -v snakemake`, i.e. the env had to be
  activated on PATH first. Worse, in a CSSwitch/macOS sandbox `$HOME` is
  virtualised to `~/.csswitch/sandbox/home`, so a `$HOME/miniforge3/...` probe
  misses the real user-level conda at `/Users/<user>/miniforge3/...`.
- **Fix**: `find_snakemake` first honours `command -v`, then probes both
  `/Users/mingleiyang` and `$HOME` across common conda roots (`miniforge3`,
  `mambaforge`, `miniconda3`, `anaconda3`) × envs (`bio-agent`, `bioagents`,
  root). The resolved `$SNAKEMAKE_BIN` is used for both the dry-run and the
  execute steps.
- **Still required**: `conda`/`mamba` must be on PATH for `--use-conda` (snakemake
  shells out to `conda info`). Prepend the conda base bin dir
  (e.g. `export PATH="/Users/mingleiyang/miniforge3/bin:$PATH"`).

## Known upstream gap (not fixed here)
- `rules/09_blast_novelty.smk` rule `blast_against_db` (and
  `rules/08_psrobot_etm.smk` rule `psrobot_etm`) require
  `results/feelnc/candidates.fa`, which **no rule in the workflow produces**
  — `gtf_to_fasta_lncrna` emits `results/feelnc/high_confidence_lncrna.fa`
  instead. With `run_blast_novelty: true` the DAG fails with
  `MissingInputException`. Keep `run_blast_novelty: false` until upstream fixes it.
