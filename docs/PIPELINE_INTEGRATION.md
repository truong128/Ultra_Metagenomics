# Scope

Ultra_Metagenomics assembles each sample once, with nf-core/mag 5.5.0 and MEGAHIT. Later stages add taxonomy, bin quality, and functional annotation. They use those contigs, proteins, and reads. They do not replace the MEGAHIT assembly.

Eukaryotic reads from a water sample are retained. The workflow does not remove them as host contamination.

| Question | Stage |
| --- | --- |
| Read quality, assembly, genes, bins | nf-core/mag 5.5.0 |
| A second read-level taxonomy | `taxprofiler` (Kaiju, MetaPhlAn, DIAMOND). Kraken2 and Bracken remain their own stages. |
| Genome-resolved metabolism | `dram2` |
| Community element cycling | `metabolic_c` |
| Read-level pathway abundance | `humann`. The nucleotide reference is built from human microbiomes, so a lake hit is encoded potential with that bias. |
| Eukaryotic contigs | `whokaryote`, plus Tiara, MetaEuk, BUSCO, and `eukcc` where those stages are enabled |
| Viral host prediction | `iphop`, after geNomad and CheckV. A prediction is not a demonstrated infection. |

## Interpretation

A sequence is encoded potential. A gene call is not evidence that the gene was active on the sampling day. An average-nucleotide-identity hit is not a GTDB species name. One sample, or one collection date, does not support a seasonal claim or a co-occurrence network.

Field metadata, blanks, chemistry, qPCR, RNA, protein, and any statistic that needs independent replicates are outside this command. A missing tool or database is recorded as BLOCKED. That status means the stage did not run.
