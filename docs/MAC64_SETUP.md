# MacBook Pro, 64 GB RAM

## A. Muc tieu / Purpose

This guide configures a resource-aware Mac mode without deleting any scientific steps from the previous 83-step registry. Core assembly, read mapping, ORF calls, two independent binners (MetaBAT2 and MaxBin2), DAS Tool refinement, CheckM2, and the optional downstream modules remain available. SemiBin2, metaSPAdes, and pplacer are retained as opt-in research/HPC paths, not started in the Mac64 default core. No GPU/CUDA required for these standard methods.

**Caveat:** `--max_memory`, Nextflow `resourceLimits`, and Docker memory allocation do **not** guarantee that every process fits 40 GB RAM. Never promise a full 45-million-pair metagenome will assemble with 64 GB until actual peaks are measured. Repeated OOM = move that stage to Linux/HPC (or change scientific design), not fabricate outputs.

## B. Installation

1. macOS recent enough for Nextflow + Java + Docker; host disk on APFS/SSD recommended. Working directory should not be mounted on an exFAT drive if Linux container permission/symlink semantics cause issues.
2. Install Java supported by the pinned Nextflow release, Nextflow, Docker Desktop (or suitable runtime), Python 3.10+.
3. Docker Desktop: Settings > Resources > Memory: **48 GB** (target; reserve ~16 GB for macOS and other work). Set CPU near 8, depending on your physical core count. Close other memory-intensive apps.
4. Apple Silicon (M-series): some `linux/amd64` bioinformatics images have no native `arm64` counterpart. Rosetta or QEMU emulation may be necessary. Check each image in an nf-core/mag test profile before a production run.
5. Edit `config/samples.tsv` to absolute FASTQ paths and `config/project.json` to the database folders. For this profile set `resource_mode` to `mac64`, with `threads=8`, `--max_cpus 8`, `--max_memory 40.GB`, and `profiles/mac64.config`. SemiBin2 is not launched by default.
6. Reference databases and sample FASTQ files are not in the repository.

## C. Quick start

```bash
cd Ultra_Metagenomics
python3 metax.py plan
python3 metax.py doctor
python3 metax.py validate --full
python3 -m unittest discover -s tests -v
python3 metax.py core --dry-run
python3 metax.py core            # launches actual core only after checking inputs/software
python3 metax.py bridge          # requires real nf-core/mag outputs
python3 metax.py audit
python3 metax.py mag-review      # creates config/approved_mags.tsv, never overwrites
```

**IMPORTANT:** `metax.py core --dry-run` only prints the Nextflow command. It is not a Nextflow pipeline test and cannot validate container compatibility or runtime RAM.

## D. MAC taxonomy from MAGs: conservative, staged

1. Review per-MAG **CheckM2** completeness and contamination metrics, confirm chimera/contamination checks and domain assignment. Mark `approved=true` in `config/approved_mags.tsv` ONLY after review and fill all metrics. The initial review sheet deliberately has `approved=false` for every MAG.
2. Strict default entry gate: completeness >=50%, contamination <=10%, GUNC `pass` or `not_run`. If GUNC not run, record incomplete validation; this is not a publication-grade purity guarantee. More conservative quality gates can be applied manually.
3. Download a skani reference database of *representative genome sketches* or build one using a version-matched representative FASTA release. This is a separate source from the GTDB-Tk package. skani does not automatically know GTDB taxonomy or its per-species ANI circumscription radii.
4. Use current skani (recommend v0.4.x when matching its v0.4 database format). Example if you already have a list of reference FASTAs (one path per line):

```bash
python3 metax.py skani-index \
  --reference-list /ABS/PATH/gtdb_reference_files.txt \
  --out /ABS/PATH/skani_gtdb_r232
```

5. Set `reference_databases.skani_representatives` to the completed DB directory; then:

```bash
python3 metax.py taxonomy-ani --dry-run
python3 metax.py taxonomy-ani
```

6. Inspect `results/taxonomy_ani/skani_search.tsv` and `results/taxonomy_ani/ani_candidates.tsv`. For every approved MAG, the latter includes ANI, AF(reference), AF(query), QC metrics and one of:
   - `ANI_CANDIDATE_REVIEW`: ANI >=95%, AF both directions >=50%; this is **only a candidate match**, not an official GTDB species call.
   - `NO_CONFIDENT_ANI_ASSIGNMENT`: there is a match but conservative thresholds not passed.
   - `NO_DATABASE_HIT`: no reference hit; does NOT imply the MAG is novel.
7. Keep unresolved MAGs in downstream analyses, at the supported higher taxonomic rank, and flag for more rigorous tree/marker-based validation. Species calls under GTDB depend on specific reference species ANI circumscription and other evidence. Do not blindly assign a species because ANI >=95%.

### Optional high-memory GTDB-Tk 2.8.0 branch

GTDB-Tk's modern `classify_wf` performs its own skani ANI screen first, but **MAGs not assigned there can still require pplacer**. Its current bacterial reference tree has a large RAM requirement, so the wrapper never silently launches it under `mac64`.

When transferred to a Linux HPC machine with adequate RAM, set `resource_mode` to `hpc` in `config/project.json`, set `reference_databases.gtdbtk` to a matching **R232** database and update all paths. Explicitly run:

```bash
python3 metax.py taxonomy-hpc                        # command preview only
python3 metax.py taxonomy-hpc --run --acknowledge-high-memory
```

The HPC command uses `gtdbtk classify_wf`, `--batchfile`, `--pplacer_cpus 1` and `--scratch_dir`, never `--full_tree`. **Do not** run the full branch on a Mac merely because an earlier skani pass looked fast. GTDB-Tk can be run in its own supported environment with matching version and DB.

## E. MAC mode computational changes from v2.0

- `resource_mode=mac64` is default, controlled by `resources` in `config/project.json`.
- `metax.py core` generates `--skip_gtdbtk`; setting `core_enable_gtdbtk=true` on Mac is an explicit error.
- Uses `--skip_spades`, `--skip_semibin`, `--skip_concoct`, `--skip_comebin`, `--skip_metabinner`. MetaBAT2 and MaxBin2 remain for DAS Tool. These other binners are not deleted: they may be re-enabled in a review/benchmark/HPC mode.
- Nextflow core cap: `--max_cpus 8 --max_memory 40.GB -c profiles/mac64.config`; same profile caps each process resource request. A local executor's actual consumption can still exceed the request; monitor Activity Monitor/Docker memory.
- Limits Nextflow JVM via `NXF_OPTS=-Xms1g -Xmx4g` unless already set by user.
- `module dram2`, `module fegenie`, `module metabolic_c`, `module antismash` now require an explicit MAG review sheet instead of treating DAS Tool output as automatically curated.
- `module eggnog`, `module dbcan`, `module dram2`, `module metabolic_c` may exceed 40 GB with certain DBs: run sequentially, reduce worker threads, change datasets/DB where scientifically valid, or offload to HPC. The external adapters do **not** enforce a 40 GB hard limit on their own.
- The previous detailed functional, eukaryotic, viral, C/N/S/P/Fe, contamination, multiomics and publication research steps are preserved in `docs/STEP_REGISTRY.tsv`; several remain conceptual or manually curated, NOT all executable in one command.

## F. Citations and documentation checked (2026-11-01)

- nf-core/mag v5.5.0 usage: https://nf-co.re/mag/5.5.0/docs/usage
- nf-core/mag v5.5.0 output: https://nf-co.re/mag/5.5.0/docs/output
- nf-core/mag parameters: https://nf-co.re/mag/5.5.0/parameters
- GTDB-Tk 2.8.0 classify_wf: https://ecogenomics.github.io/GTDBTk/commands/classify_wf.html
- GTDB-Tk resource requirements: https://ecogenomics.github.io/GTDBTk/installing/
- skani repo (current v0.4 sketch format compatibility): https://github.com/bluenote-1577/skani
- skani reference indexing: https://github.com/bluenote-1577/skani/wiki/Tutorial:-setting-up-the-GTDB-genome-database-to-search-against
- Nextflow resource limits: https://docs.seqera.io/nextflow/reference/process/directives/resource-limits
- ARM architecture note: https://nf-co.re/blog/2026/arm-pipelines
