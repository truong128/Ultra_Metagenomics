# Retained previous stages: 24-step crosswalk

Everything in the previous analytical plan is mapped; nothing disappears merely because different tools overlap.

| Old | Previous analytical scope | New implementation/work package |
|---|---|---|
| 00 | Sampling, environmental metadata and controls | Registry 00-09 |
| 01 | FASTQ intake and integrity | 10-11 (local validator) |
| 02 | FastQC, fastp, MultiQC | 12-13 (core) |
| 03 | PhiX, laboratory and host-specific contamination | 14-17 (core and curated controls) |
| 04 | Read-level taxonomic classification | 18-20 (Kraken2; Bracken/Kaiju optional) |
| 05 | SSU ribosomal 16S/18S phyloFlash | 21-22 (module and tree QC) |
| 06 | Assembly and coassembly strategy | 23-26 (production vs bench modes) |
| 07 | QUAST, Bowtie2, mapping and ALE | 27-30 (core) |
| 08 | CAT/Tiara contig taxonomy | 60-62 + core optional CAT |
| 09 | Prodigal, MetaEuk, tRNA/rRNA | 31-33 (core and optional) |
| 10 | Nonredundant gene catalogue and abundance | 34-35, 71 |
| 11 | eggNOG, KOfam, Pfam, UniRef, EC, CAZyme | 36-41 |
| 12 | Pathway reconstruction | 48-56 |
| 13 | MAG ensemble, refining, quality, taxonomy | 42-47 |
| 14 | DRAM2, METABOLIC-C, MAG function | 48-50 |
| 15 | Eukaryotic validation and euk-MAG | 60-62 |
| 16 | Virus, phage, provirus, host linkage | 63-66 |
| 17 | CAZyme and substrate use | 39-40, 52 |
| 18 | C/N/S/P/salinity/iron/metal/stress/motility | 50-59 |
| 19 | AMR, BGC, CRISPR, virulence | 67-70 |
| 20 | Read-/MAG-/gene-level quantification | 71 |
| 21 | Diversity and multivariate/interaction statistics | 72-75 |
| 22 | RNA, protein, wet-lab metabolic validation | 76 |
| 23 | Figures, methods, metadata, deposition, QA | 77-82 |

A crosswalk is not proof each analytical result can be obtained from a single sample or that every listed extension is fully automated. Consult `STEP_REGISTRY.tsv` statuses.
