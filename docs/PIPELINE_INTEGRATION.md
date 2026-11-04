# How the external pipelines were merged

The shared design chat compared six whole workflows and a set of specialist tools. Running all of them on the same reads would assemble, call genes, and count abundance more than once, and those counts cannot be added together. v3.1 keeps one production engine and turns the other workflows into either a non-overlapping evidence stage or an explicit benchmark.

| Source | What it is good at | What v3.1 does with it |
|---|---|---|
| [nf-core/mag 5.5.0](https://nf-co.re/mag/5.5.0/) | QC, MEGAHIT, mapping, Prodigal, ensemble binning, CheckM2, BUSCO, Tiara, geNomad, MetaEuk, GUNC | Production engine. Full binner set on linux64. |
| [nf-core/taxprofiler](https://nf-co.re/taxprofiler) | Kraken2, Bracken, Kaiju, MetaPhlAn, DIAMOND from reads | Stage `taxprofiler`, after a database CSV is set. Kraken2/Bracken remain the existing read-taxonomy stages so they are not launched a second time inside taxprofiler unless the CSV asks for them. |
| [SqueezeMeta](https://github.com/jtamames/SqueezeMeta) | One-command assembly, genes, taxonomy, function, MAGs | Stage `benchmark_squeezemeta`. **Held** until `benchmark_engines` contains `squeezemeta`. Output stays in its own directory. |
| [DRAM / DRAM2](https://github.com/WrightonLabCSU/DRAM) | Genome-resolved metabolism | Existing `dram2` stage, unchanged role. |
| [METABOLIC](https://github.com/AnantharamanLab/METABOLIC) | Community C/N/S contributions | Existing `metabolic_c` stage, unchanged role. |
| [AAFC metagenomics-snakemake](https://github.com/AAFC-Bioinfo-AAC/metagenomics-snakemake) | QC, taxonomy, KEGG, MAGs, AMR, dbCAN for host-associated samples | Scope already covered by nf-core/mag, eggNOG/KOfam, AMRFinderPlus and dbCAN. Host depletion is not imported: water samples keep eukaryotic DNA. |
| [HUMAnN](https://github.com/biobakery/humann) | Read-level pathway abundance | Stage `humann`. Results are encoded potential. Its nucleotide reference is biased toward human microbiomes, so it does not replace eggNOG or METABOLIC for mangrove water. |
| [metagenome-atlas](https://github.com/metagenome-atlas/atlas) | Alternative Snakemake MAG engine | Stage `benchmark_atlas`, held the same way as SqueezeMeta. |
| Whokaryote | Eukaryote vs prokaryote contigs from gene structure | Stage `whokaryote`, model S, so it does not re-run Tiara. |
| EukCC, BUSCO, MetaEuk | Eukaryotic genes and MAG quality | BUSCO inside the core. MetaEuk when its database exists. EukCC only on bins already called eukaryotic. |
| iPHoP | Viral host prediction | Stage `iphop` after CheckV/geNomad. A prediction is not a demonstrated infection. |
| dRep, minced, MMseqs2, KOfamScan, Nonpareil, HMMER, CoverM, skani, GTDB-Tk, FeGenie, antiSMASH, CheckV | Specialist steps already listed in the registry | Each is its own stage. None repeats MEGAHIT or Prodigal. |

## Held on purpose

Set this only when a second assembler is a benchmark, and do not mix its contigs into the production audit:

```json
"benchmark_engines": ["squeezemeta", "atlas"]
```

metaSPAdes stays off unless `core_benchmark_spades` is true. That flag removes `--skip_spades` from the same nf-core run. It is a comparison, not a second biological result.

## Still not claimed

Field metadata, blanks, environmental chemistry, qPCR, RNA, protein, and any statistic that needs independent replicates are still registry steps the FASTQ command cannot finish. A KO or HMM hit is encoded potential, not measured activity. ANI candidates are not GTDB species names. One water sample cannot support a seasonal or network claim.
