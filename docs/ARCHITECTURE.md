# Architecture

## Production assembly

One production engine: **nf-core/mag v5.5.0**, Nextflow/Docker, per-sample Illumina paired-end assembly. The core performs QC, fastp, PhiX removal, MEGAHIT assembly, QUAST, read mapping, gene prediction, binning, bin QC, and optionally GUNC, Tiara, geNomad and MetaEuk. Production contigs are the MEGAHIT contigs. Downstream stages use those contigs and the Prodigal proteins from the same run.

Domain-specific evidence is added after the core, using the **same contigs and Prodigal proteins**, so the reads are not assembled or called a second time. Strict `assets.json` bridging provides explicit sample identifiers, assembly FASTA, protein FASTA, Prodigal GFF and DAS-refined candidate MAG FASTAs. Its source paths correspond to the official output specification for nf-core/mag 5.5.0: `Assembly/MEGAHIT/<sample>.contigs.fa.gz`, `Annotation/Prodigal/MEGAHIT-<sample>.faa.gz` and `.gff.gz`. It does not guess filenames from unknown source structures. It rejects missing files. It refuses coassembly where the one-to-one sample bridge would be inaccurate. When DAS Tool yields no bins, the manifest warns but does not misreport zero MAGs.

**Production execution graph:**

```
metadata, controls, raw FASTQs
       |
validate_inputs -> nf-core/mag 5.5.0 [FastQC, fastp, PhiX, MEGAHIT, QUAST, ALE, mapping, Prodigal, binners, DASTool, CheckM2, Tiara bins]
       |                                       |
       +--> MAG/bin_summary/prodigal outputs <-+
       |                                       |
       +--> strict bridge: assets.json --------+
                 |                 |               |
        SSU read taxonomy       ORFs              MAGs
        phyloFlash             eggNOG            DRAM2
        Kraken2 (+ Bracken)    dbCAN             FeGenie
        protein Kaiju         KOfam             antiSMASH
        PR2/phylogenetics     AMRFinderPlus      GTDB (core optional)
                \              |                 /
                 ---> curated metabolic ledger <---
                             |
                     METABOLIC-C (community)
                             |
                       C/N/S/P, iron, salinity, stress
                             |
                environmental statistics (only with replication)
                             |
                    figures, biological validation,
                    MIxS/MIMAG/MIUViG, FAIR deposition
```

## De-duplication table

| Function | Default | Alternative | Decision |
|---|---|---|---|
| Quality + trimming | nf-core/mag fastp + FastQC | standalone fastp | Standalone QC is NOT run in production |
| Assembly | nf-core/mag MEGAHIT | a second assembler | Production tables use the MEGAHIT contigs |
| Gene calls | nf-core/mag Prodigal | a second gene caller | Reuse the core proteins |
| Read mapping | nf-core/mag Bowtie2 | custom mapper | No repeated mapping unless extra scientific question |
| Binning | MetaBAT2 + MaxBin2 (Mac64); optional SemiBin2 research/HPC | others | Different algorithms genuinely add genome recovery; DAS Tool outputs curated bins |
| MAG QC | CheckM2 + optional GUNC | CheckM1 | No routine CheckM1/2 duplicates |
| Taxonomy | Kraken2 read-level, GTDB-Tk MAG-level | Kaiju read-level | Different biological unit/evidence, not redundant |
| Euk | Tiara on contigs + SSU + MetaEuk on bins | PR2/18S phylogeny | Complementary rank and gene model |
| Viral | geNomad in nf-core/mag | standalone geNomad | Never rerun if core enabled; CheckV adds genome quality |
| Gene functions | eggNOG-mapper ORF | DRAM2 (MAG-level) | Different resolution, cache gene IDs and compare consistency |
| CAZyme | dbCAN dedicated | DRAM2 dbCAN | Use dedicated only for CGC/substrate or gene catalogue; report overlap |
| C/N/S pathways | curated KO/HMM + METABOLIC-C | pathway labels from eggNOG alone | Avoid generic gene hits interpreted as complete pathways |

## Environmental-specific hazards

- A single water sample cannot yield unbiased biodiversity richness or seasonal dynamics; contig count != species count.
- Environmental eukaryotic DNA is signal, not automatically contaminant. Do not run blanket eukaryote host depletion.
- Read taxonomy can use raw reads in the provided standalone module but this must be reported as **raw-read profiling**, not identical denominator to cleaned assembled reads. To profile trimmed reads use the `--save_clipped_reads` outputs and a manually checked read manifest. nf-core/mag also removes PhiX downstream; pay attention to differences.
- `nifH`, `amoA`, `pmoA`, `dsrAB`, `sox`, `nosZ`, `rbcL`, `cyc2` have paralogs/context/phylogenetic issues. Validate marker families and directionality with gene trees, neighborhood and MAG support.
- Full N/S/C pathways assigned to the entire community must not be assigned to individual organisms unless a MAG supports their completeness.
- For viral sequences use validated viral reference genomes, provirus context, CheckV; a candidate vOTU != a named virus species.
- Negative controls allow contamination assessment, not automatic removal of taxa without evaluation.
- For comparisons: compositional methods, batch covariates, FDR, effect sizes, pre-registration, independent samples. Network connections are hypotheses.
- RNA, proteomics, metabolomics and environmental chemistry are distinct experimental evidence. DNA alone indicates **encoded potential**, not in-situ metabolic flux.

## Robustness: error barriers

- Input parser checks paired headers, FASTQ structure, read length vs quality length and real file existence.
- Explicit tool selection, version tag and false-positive sensitive thresholds; optional DBs never silently substituted.
- Bridge does not infer absent MAG quality, species, taxon or read percentages from contig counts.
- Built-in audit calculates assembly/protein FASTA metrics; no made-up biology.
- Downstream scripts abort on missing mandatory DB/program and never silently mark completion.
- Publication gate holds for non-replicated designs and unvalidated marker claims.

## Additional research-grade extensions not yet fully automated here

Planned rather than runnable-with-one-command: marker phylogenetic trees, reaction-specific ecological pathway reconstruction, model-based contamination controls, genotype-resolved euk-MAG binning, spike-in absolute quantification, vOTU coverage/host-prediction, HGT, strain microdiversity, protein structures, redox-aware model, pathway flux measurements, RNA/protein evidence integration, spatiotemporal statistics and deposition upload. They remain documented in `STEP_REGISTRY.tsv`; this prevents claiming that every conceptual step is implemented.

`--postbinning_input refined_bins_only` means CheckM2/BUSCO evaluate only DAS Tool refined bins in the primary run. This avoids scoring redundant raw bins several times and is consistent with the bridge, but requires a fallback for samples with no passing refined bins. If a refined bin is absent, inspect original binner candidates and repeat under a versioned fallback config; never silently invent a zero-MAG outcome.

## v2.1 local ANI-first taxonomy and QC gate

Default Mac64 mode disables GTDB-Tk inside nf-core/mag. It is not enough to set a single `--cpus` option, since its `pplacer` phase can require much more than 64 GB RAM. Approved MAG FASTAs (quality reviewed using CheckM2, optionally GUNC) feed standalone skani against **separately prepared** representative-genome sketch DB. `ani_candidates.tsv` reports ANI and AF and does not call GTDB species. Unassigned remains unassigned, with a high-RAM GTDB-Tk v2.8.0 branch explicitly available only for HPC.

Mac64 Nextflow `-c profiles/mac64.config` + `--max_memory 40.GB --max_cpus 8` cap per-task *requested* resources; this does not control physical peak consumption of every local task or every downstream module. The Mac package keeps all prior functional branches. `docs/MAC64_SETUP.md` is authoritative for this profile.
