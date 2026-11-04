# MangroveMetaX Ultra changelog

## v3.1 Linux workstation (2026-11-02)

- Default `resource_mode` is `linux64`: 24 CPUs, 48 GB per container, one heavy task at a time, NVIDIA runtime visible to containers.
- Full binning ensemble (MetaBAT2, MaxBin2, SemiBin2, CONCOCT, COMEBin, MetaBinner) plus BUSCO. COMEBin and MetaBinner may fail without discarding the other binners.
- geNomad, MetaEuk, GUNC and CAT attach only when their database path exists, so assembly is not blocked by a missing optional database.
- GTDB-Tk pplacer is a separate `taxonomy_gtdb` stage (1 pplacer CPU, disk scratch, no full tree). It stays off until RAM+swap is at least 96 GB.
- Added non-overlapping stages for taxprofiler, Whokaryote, Nonpareil, KOfamScan, MMseqs2, CoverM, dRep, EukCC, iPHoP, minced, HUMAnN and a biogeochemical HMM screen.
- SqueezeMeta and metagenome-atlas are held benchmark engines. The AAFC host-associated workflow is not run, and host depletion is not imported.
- Mac64 behavior is unchanged and now lives in `config/project_mac64.json`.

## v3.0 AutoRun (2026-11-01)

- Added `metax.py run-all` with topological ordering, resume, `--dry-run`, `--force`, strict completion semantics, opt-in `--allow-partial`.
- Added KO annotation-count summary, measured-FASTA SVGs, per-stage stdout/stderr logs and status.json, report.json, report.html with explicit data/tool/DB blockers and incomplete research steps.
- Added automatic conservative MAG candidate QC from nf-core/mag CheckM2 summary (exact IDs, hard stop on mismatch, no false official taxonomy).
- Automated Bracken from Kraken2 reports, adapted CheckV to multiple viral FASTAs and decompression of METABOLIC-C input FASTQ.
- Preserved Mac64 guard rails, explicit high-memory GTDB-Tk HPC workflow, scientific boundary warnings.
- All previous research steps retained in the 91-step registry; some remain non-computational/manual rather than falsely reported as automated.
- Local unit tests added to check automated QC, fail-closed behavior, scheduler preflight and truthful incomplete reports.

## v2.1

Mac64 mode and skani ANI candidate screening; v2.1 README replaced by v3.0 instructions.
