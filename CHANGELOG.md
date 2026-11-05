# Changelog

## 2026-11-05

- Public documentation describes the production workflow. MEGAHIT, through nf-core/mag 5.5.0, is the production assembler.
- HUMAnN 3.9 uses MetaPhlAn 4.1.1 and the June 2023 markers. The January 2026 MetaPhlAn database stays with read-level taxonomy.

## 2026-11-02

- Default `resource_mode` is `linux64`: 24 CPUs, 48 GB per container, one heavy task at a time.
- Binning uses MetaBAT2, MaxBin2, SemiBin2, CONCOCT, COMEBin, and MetaBinner, then DAS Tool, CheckM2, and BUSCO. COMEBin and MetaBinner may fail without discarding the other bins.
- geNomad, MetaEuk, and GUNC run when their database path is set.
- GTDB-Tk is a separate `taxonomy_gtdb` stage, outside the 48 GB Nextflow cap, with one pplacer CPU and no full tree. It starts when RAM plus swap is at least 150 GB.
- Additional stages cover a second taxonomy, eukaryotic contigs, pathway abundance, virus-host prediction, and marker screens. Each stage uses the MEGAHIT contigs or the original reads. None of them assembles the reads again.
- The Mac profile is unchanged and lives in `config/project_mac64.json`.

## 2026-11-01

- Added `metax.py run-all` with resume, `--dry-run`, `--force`, and an opt-in `--allow-partial`.
- Added an execution report, per-stage logs, a KO count summary, and FASTA summary figures.
- MAG screening uses the CheckM2 summary at 50% completeness and 10% contamination. A failed match stops the gate. The gate does not assign a species name.
- Bracken runs from the Kraken2 report. CheckV accepts the viral FASTA from geNomad.
- Steps that are not computed from FASTQ stay in `docs/STEP_REGISTRY.tsv` and are not reported as finished.
