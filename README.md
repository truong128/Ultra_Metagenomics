# Ultra_Metagenomics

Shotgun metagenome analysis for water samples. The production assembly is nf-core/mag 5.5.0 with MEGAHIT, one sample at a time. Later stages add taxonomy, genome bins, and functional annotation. They do not replace the MEGAHIT contigs.

The repository contains the workflow code and a sample sheet. It does not contain sequence reads, reference databases, or results. Those files are large, and the reads are already public at ENA.

The workflow is written for Ubuntu. A MacBook Pro can run a defined subset. See Hardware.

## Worked example

The sample sheet is two Lake Motosu libraries from the freshwater lake study PRJEB35640 (ERP118719). Both were collected on 7 September 2015. One date and two depths support a depth comparison. They do not support a seasonal claim or a co-occurrence network.

| Sample | Depth | ENA run | BioSample | Alias |
| --- | --- | --- | --- | --- |
| Motosu_5m | 5 m | ERR14787630 | SAMEA117850005 | Mot-5 |
| Motosu_100m | 100 m | ERR14787631 | SAMEA117850006 | Mot-100 |

The same study contains other European and Asian lakes. They are not in `config/samples.tsv`. Add a row only after the matching FASTQ files are on disk.

Download the two Motosu pairs into the paths named in the sample sheet. Each file is about 1.1 GB.

```bash
mkdir -p data/ERP118719
curl -L -o data/ERP118719/ERR14787630_1.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/030/ERR14787630/ERR14787630_1.fastq.gz
curl -L -o data/ERP118719/ERR14787630_2.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/030/ERR14787630/ERR14787630_2.fastq.gz
curl -L -o data/ERP118719/ERR14787631_1.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/031/ERR14787631/ERR14787631_1.fastq.gz
curl -L -o data/ERP118719/ERR14787631_2.fastq.gz \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR147/031/ERR14787631/ERR14787631_2.fastq.gz
```

`config/accessions.tsv` repeats the accessions and the FTP paths. FASTQ files stay out of git (`data/` is ignored).

## Two ways to run

Install Nextflow, Docker, and the databases you intend to use. Copy a profile and point it at those databases:

```bash
cp config/project.linux64.json config/project.json
# edit reference_databases in config/project.json
# databases live outside the repo; or set METAX_DB_ROOT
python3 metax.py doctor
python3 metax.py run-all --dry-run
python3 metax.py run-all
```

`run-all` is sequential and resumable. Run the same command again after a stop. A missing tool or database is recorded as BLOCKED. That status means the stage did not run. It is not evidence that a gene, taxon, or virus is absent. `--allow-partial` writes the report anyway and is not a finished analysis.

For a Mac, copy `config/project_mac64.json` instead and stay inside the limits in that file: 8 threads and 40 GB per process.

Empty strings in `reference_databases` turn the matching stages off. Fill a path only for a stage you want. The groups below say which paths a water-sample study should fill.

### Research track

Use this track for a lake or other water metagenome when the questions are community composition, genome bins, and the main encoded functions. It is enough for those questions. Stages outside this list can be left empty.

Run the groups in this order. Within a group, the tools answer the same kind of question and can share one database choice where the table says so.

**1. Reads and contigs.** Required.

| Stage | What it answers |
| --- | --- |
| `validate` | The pairs are intact FASTQ. |
| `core` QC | fastp trims the reads. FastQC records quality. These reads are about 150 bp, which is the Bracken read length in the config. |
| MEGAHIT | Assembles each sample. metaSPAdes stays off so a second assembler cannot enter the contig tables. |
| Mapping and QUAST | Places reads back on the contigs and records assembly size. |
| Prodigal | Predicts prokaryotic genes on the contigs. |

**2. Who is in the water.** Required for a community description.

| Stage | What it answers | Database to install |
| --- | --- | --- |
| `taxonomy` | Read-level taxonomy, independent of the assembly. | Kraken2 PlusPFP 16 GB index. About 16 GB on disk. |
| `bracken` | Corrects that taxonomy toward read abundance. | Same Kraken index. |
| `ssu` | Reconstructs SSU rRNA, which is the marker most water papers still report beside shotgun taxonomy. | SILVA release used by phyloFlash. About 16 GB. |

**3. Genome bins, if the paper uses MAGs.** Required for that paper. Skip the whole group if the study stays at reads and contigs.

| Stage | What it answers |
| --- | --- |
| MetaBAT2, MaxBin2, CONCOCT | Three bin sets from the same contigs. |
| DAS Tool | One refined bin set. |
| CheckM2 | Completeness and contamination. |
| `mag_auto` | Keeps bins at or above 50% completeness and at or below 10% contamination. This is an automatic gate. It is not a manually curated species list. |
| `coverm` | Abundance of those bins from the existing BAM files. |
| `drep` | Collapses the same population recovered from both depths so it is not counted twice. Worth running once there is more than one sample. |
| GUNC | Screens refined bins for chimeras. Include it when the bins will be published. |

SemiBin2 can be added on Ubuntu. It is off in the Mac profile. COMEBin and MetaBinner are allowed to fail without discarding the other bins. Leave them for the full track unless you have already budgeted the extra time and memory.

**4. Encoded functions.** Required: one ortholog database, plus dbCAN for a water sample.

| Stage | What it answers | Install |
| --- | --- | --- |
| `kofam` or `eggnog` | KEGG Orthologs, or the broader eggNOG orthologs. One of the two is enough. KOfam is the smaller install (about 9 GB). eggNOG is about 48 GB. | One path. |
| `dbcan` | Carbohydrate-active enzymes. Polysaccharide degradation is a standard question in lake water, so this stage belongs in the research track even when the ortholog stage is already present. | dbCAN, about 8 GB. |
| `functional_summary` | Counts ortholog genes per sample from the eggNOG table. It runs after `eggnog`. The counts are gene presence in the assembly, not activity and not transcript abundance. | None beyond eggNOG. |

**5. Viruses, when the viral fraction is part of the question.** Recommended for lake water. Optional when the paper is prokaryotic MAGs only.

| Stage | What it answers | Install |
| --- | --- | --- |
| geNomad, inside `core` | Viral sequence on the contigs. | nf-core/mag downloads this database when `core_enable_genomad` is true. |
| `checkv` | Completeness and contamination of those viral sequences. | CheckV database, about 11 GB. |

**6. Species names on the bins.** Required only when a species name will appear in the paper.

| Stage | What it answers |
| --- | --- |
| `taxonomy_gtdb` | GTDB-Tk `classify_wf` on the bins that passed the CheckM2 gate. This is the species call. |
| `taxonomy_ani` | skani against GTDB representative sketches. A hit is an ANI candidate. It is not a GTDB species assignment. Use it when GTDB-Tk will not run, including every Mac run. |

GTDB-Tk release 232 is about 94 GB after extraction. The stage needs about 150 GB of RAM plus swap and does not run inside the 48 GB Nextflow cap. On 64 GB of RAM it needs a swap file. See Hardware.

**7. Bookkeeping.** Always produced.

`bridge` links contigs, proteins, and bins for the later stages. `audit` records assembly and protein metrics. `figures` draws three FASTA summaries. `report` writes `results/execution/report.html`.

### Full track

Everything in the research track, plus the stages below. Each one answers an extra question. A water study is complete without them. `run-all` runs a stage only when its program is on `PATH` and its database path is set. Leave the path empty to skip it.

**More taxonomy, still not required.**

| Stage | Extra question | Disk |
| --- | --- | --- |
| `taxprofiler` | Kaiju, MetaPhlAn, and DIAMOND/SwissProt as a second read-level taxonomy. | Kaiju reference plus fungi is about 111 GB. MetaPhlAn is about 85 GB. |
| `taxonomy_full` | Kraken2 on the uncapped PlusPFP index, memory-mapped, after every other stage. The 16 GB index remains the one used by `taxonomy`. | The compressed archive is 200 GB. The extracted index needs about 270 GB. Keep the 16 GB index either way. |

**Eukaryotes in the water column.** Lake water contains protists and fungi. Add this group when they are part of the question. Do not strip eukaryotic reads before assembly if you might add it later.

| Stage | Extra question |
| --- | --- |
| `euk_contig` (Tiara) and `whokaryote` | Which contigs look eukaryotic from sequence composition or gene structure. |
| MetaEuk, inside `core` | Eukaryotic gene prediction against UniRef50. |
| BUSCO, inside `core` | Marker completeness, with automatic lineage selection. |
| `eukcc` | Completeness of eukaryotic bins. The database is about 18 GB. |

**Viruses, beyond geNomad and CheckV.**

| Stage | Extra question | Disk |
| --- | --- | --- |
| `iphop` | Predicted host of a viral sequence. | About 294 GB. |
| `crispr` | CRISPR arrays on the QC-threshold MAGs. | Tool only. |

**Narrow metabolic and product questions.**

| Stage | Extra question | Disk |
| --- | --- | --- |
| `markers_hmm` | A short HMM screen for carbon, nitrogen, and sulfur marker genes. Useful in water. Redundant once KOfam or eggNOG already covers those pathways. | Small HMM file. |
| `dram2` | A second metabolic summary of the bins (KOfam, dbCAN, MEROPS, VOG, FeGenie HMMs, CAMPER, UniRef). KEGG is not included. The Nextflow run needs Java 17 through 26. | About 102 GB, of which UniRef90 is 96 GB. |
| `fegenie` | Iron genes. Reasonable for a deep, stratified lake if iron is a hypothesis. | Tool data only. |
| `metabolic_c` | Community metabolism with METABOLIC-C. | Small, but it unpacks reads and uses a lot of scratch. |
| `antismash` | Biosynthetic gene clusters on the bins. | antiSMASH databases, separate from the table above. |
| `amrfinder` | Antimicrobial-resistance gene candidates. A specific question, not a default lake result. | About 0.2 GB. |
| `humann` | Read-level pathway abundance. HUMAnN 3.9 runs with MetaPhlAn 4.1.1 and the June 2023 markers `mpa_vJun23_CHOCOPhlAnSGB_202403`. The January 2026 MetaPhlAn database stays with `taxprofiler` and is not a valid HUMAnN input. The nucleotide reference is built from human microbiomes, so lake hits stay encoded potential with that bias. | About 50 GB, plus about 25 GB for the June 2023 markers. |
| `catalogue` | A non-redundant protein catalogue of the predicted genes. | Scratch for MMseqs. |
| `nonpareil` | How redundant the read set is, and a rough sequencing-effort curve. | Tool only. |

## Hardware

**Ubuntu, research track.** 64 GB RAM is enough for QC, MEGAHIT, mapping, MetaBAT2, MaxBin2, CONCOCT, DAS Tool, CheckM2, Kraken 16 GB, Bracken, phyloFlash, KOfam or eggNOG, dbCAN, and CoverM. Keep at least 1 TB free. With KOfam, SILVA, dbCAN, CheckV, and the 16 GB Kraken index, the research databases are under 120 GB. Assembly and mapping need on the order of 100 GB of scratch for libraries of this size. For two libraries like the worked example, a full run can use about 300 GB for Nextflow work and published results.

**Ubuntu, full track.** 64 GB RAM is the minimum: a 48 GB ceiling on each Nextflow container, and one heavy container at a time. GTDB-Tk does not fit in 64 GB. Provide 128 GB RAM, or 64 GB RAM plus a swap file of at least 128 GB. The profile starts that stage only when RAM plus swap is at least 150 GB. Swap is slower than RAM. The uncapped Kraken index is memory-mapped, so it does not need 270 GB of RAM. Run it after the other stages. Keep 4 TB free before installing the full database set.

| Database | Size | Track |
| --- | --- | --- |
| iPHoP | 294 GB | Full |
| Kraken, full index | 270 GB extracted | Full. Keep the 16 GB index either way. |
| GTDB release 232 | 94 GB | Full, for species calls |
| Kaiju | 111 GB | Full |
| DRAM2 | 102 GB | Full |
| MetaPhlAn | 85 GB | Full |
| HUMAnN | 50 GB | Full, plus about 25 GB for the June 2023 markers |
| eggNOG | 48 GB | Research, when KOfam is not used |
| EukCC | 18 GB | Full |
| SILVA / phyloFlash | 16 GB | Research |
| Kraken PlusPFP-16 | 16 GB | Research |
| CheckV | 11 GB | Research, when viruses are included |
| KOfam | 9 GB | Research |
| dbCAN | 8 GB | Research |

**MacBook Pro.** Use `resource_mode=mac64`. The profile caps the run at 8 threads and 40 GB per process, skips GTDB-Tk, skips SemiBin unless `enable_semibin_on_mac` is set, and skips metaSPAdes. The shipped Mac config also leaves geNomad, MetaEuk, GUNC, and BUSCO off. A MacBook Pro with 64 GB of unified memory can run the research stages that stay inside 40 GB: QC, MEGAHIT, MetaBAT2, MaxBin2, DAS Tool, CheckM2, Kraken 16 GB, Bracken, phyloFlash, KOfam or eggNOG, dbCAN, CoverM, dRep, and Nonpareil. Do not plan GTDB-Tk, the 270 GB Kraken index, iPHoP, DRAM2 with UniRef90, or COMEBin on a laptop. Apple silicon may be running Linux amd64 containers through Rosetta, which is slower. Put databases and `work/` on an external SSD. Keep at least 1 TB free for the research track. A 16 GB or 32 GB MacBook is outside this profile.

## Reading the output

DNA sequence is encoded potential. It does not show which genes were active on the sampling day.

A skani hit at the thresholds in the config (95% ANI, 50% aligned fraction) is a candidate. The species label used in text should come from `taxonomy_gtdb` when that stage ran.

CheckM2 at 50% / 10% is the automatic bin gate. Bins outside that gate remain on disk. They are excluded from stages that ask for the QC-threshold set.

One blocked stage does not zero the taxon or pathway it would have tested. Read `results/execution/report.html` before treating a table as complete.

## Layout

```text
config/samples.tsv              sample sheet; FASTQ paths are local, files are not in git
config/accessions.tsv           ENA accessions for the worked example
config/project.linux64.json     Ubuntu profile; copy to config/project.json
config/project_mac64.json       MacBook Pro profile
profiles/linux64.config         Nextflow ceiling: 24 CPU, 48 GB, one heavy task
profiles/mac64.config           Nextflow ceiling: 8 CPU, 40 GB
docs/LINUX64_SETUP.md           swap and the GTDB-Tk stage
docs/DATABASE_SETUP.md          where each database path is documented
```

`config/project.json` is gitignored after you copy and fill it, so local paths are not pushed.

## Citation

Truong Nguyen, CUNY School of Medicine.
