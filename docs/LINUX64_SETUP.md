# Linux profile (64 GB RAM)

`config/project.linux64.json` sets `resource_mode=linux64`. Copy it to `config/project.json` and fill the database paths before a run. The profile was exercised on Ubuntu with:

- 16 cores / 32 threads
- 64 GB RAM
- an NVIDIA GPU exposed to Docker with `--gpus all`
- one heavy task at a time, 24 CPU and a 48 GB per-container ceiling

A GPU is not required for assembly, binning, or the annotation stages. The Mac profile is `config/project_mac64.json` and `profiles/mac64.config`. Disk and memory figures for a two-sample lake run are in the repository `README.md`.

## What is on now

Production engine remains **nf-core/mag 5.5.0**. MEGAHIT is still the only production assembler. The core command now also requests:

- MetaBAT2, MaxBin2, SemiBin2 (`--semibin_environment global`), CONCOCT, COMEBin and MetaBinner, then DAS Tool and CheckM2
- BUSCO with automatic lineage selection
- Tiara bin-domain classification
- geNomad, MetaEuk, GUNC and CAT only when that database path exists

COMEBin and MetaBinner are allowed to fail without discarding the other binners. They are the two most likely to exceed 64 GB.

GTDB-Tk pplacer is a separate stage, `taxonomy_gtdb`. It uses one pplacer CPU and a disk scratch directory. It does not use `--full_tree`. It does not run inside the 48 GB Nextflow cap, because that cap stops pplacer. Bacterial classification wants on the order of 100–150 GB of RAM plus swap. The linux profile blocks the stage until that total is available (`gtdbtk_min_ram_plus_swap_gb`, default 150).

On a 64 GB machine, add a swap file of at least 128 GB. A 384 GB file is what was used for the Lake Motosu run described in `README.md`. pplacer then runs from swap and is slow. The stage still needs `reference_databases.gtdbtk` and `gtdbtk` on `PATH`.

```bash
sudo fallocate -l 128G /swap-metagenomics.img
sudo chmod 600 /swap-metagenomics.img
sudo mkswap /swap-metagenomics.img
sudo swapon /swap-metagenomics.img
```

Add one matching line to `/etc/fstab` if the file should survive a reboot. Do not run `mkswap` again on a file that is already active.

## Run

```bash
cd MetagenomicsX
python3 metax.py doctor
python3 metax.py run-all --dry-run
python3 metax.py run-all
```

`run-all` is sequential and resumable. A missing tool or database is BLOCKED and written to `results/execution/report.html`. It is not reported as a biological negative. `--allow-partial` accepts that incomplete report without a nonzero exit.

## Database paths

Set absolute paths in `config/project.json` → `reference_databases`. Empty values are skipped. The core assembly still runs. See `docs/DATABASE_SETUP.md` and `docs/PIPELINE_INTEGRATION.md`.
