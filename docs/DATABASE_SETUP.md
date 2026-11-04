# Database installation & provenance

Large references are NOT included. Download each via the official upstream instructions into a local Linux filesystem, not a network/cloud or exFAT working directory. The pipeline cannot generate trusted taxonomy with incomplete or irrelevant reference databases.

- Kraken2: environmental customized bacterial/archaeal/viral/eukaryotic DB, licensed; calibrate Bracken kmer distribution for actual read length if Bracken enabled. Large euk DBs may require much more RAM than 64 GB.
- SILVA/phyloFlash: versioned SSU reference (16S, 18S). Independently verify protist placements with PR2 and phylogenetic alignments.
- eggNOG-mapper: download officially supported eggNOG annotation/alignment dataset and store `reference_databases.eggnog_data`.
- dbCAN: install dbCAN3/5 compatible database, annotate proteins and separately run CGC pipeline only when GFF/genomic context matches source.
- MetaEuk: if activated in core, provide correctly indexed MMseqs2/MetaEuk DB. nf-core/mag only runs MetaEuk when enabled.
- geNomad: optional activation inside nf-core/mag, provide a valid DB. External `checkv` requires a matched separately downloaded CheckV DB.
- DRAM2: separate plugin workflow with its own databases; this is a 2026 beta workflow and should be pinned/benchmarked. `--sum_ecos` depends on current DRAM2 rules; verify available ecosystem rules before production.
- GTDB-Tk, CheckM2, GUNC: environment-specific databases can consume hundreds of GB; download/pin only when enabled and scientifically justified. GTDB only applies to Bacteria/Archaea.
- AMRFinderPlus: update and record the NCBI pathogen database and organism settings where relevant; any AMR signal from ambient DNA does not imply clinical phenotype.
- FeGenie: requires FeGenie-specific HMM and R dependencies configured as documented by its authors.
- antiSMASH: run selected MAGs after QC (metagenome gene fragments can give partial BGCs).
- METABOLIC-C: version 4.0 upstream tool with version-specific requirements, separate env; input MAG directory and reads/BAM must match what was mapped. This is NOT silently run as part of DRAM2.

Record release, acquisition date, source/licence, hash and absolute path in `config/db_manifest.tsv`. If `CHECKV_DB` etc are empty, modules refuse to run rather than making fake reference calls.

## v2.1 skani ANI index (Mac profile)

Standalone skani index requires its own sketch directory generated from the selected representative FASTA release. Set `reference_databases.skani_representatives` to it. Do **not** assume GTDB-Tk internal skani sketch paths can be passed unchanged to standalone skani. A GTDB version, genome reference release, skani release and checksums should be recorded. `skani` v0.4 has a sketch format update and should be paired with its own v0.4 index where practical. Store `GTDB_skani_representatives` in `config/db_manifest.tsv`. Low-RAM search avoids loading the `pplacer` reference tree, but construction of the original reference index may itself use extensive disk and RAM. Pre-indexed data may be preferable on a 64 GB machine.
