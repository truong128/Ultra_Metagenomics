# Publication gates and acceptance criteria

| Gate | Automated portion | Scientific decision (cannot be automated blindly) |
|---|---|---|
| G0 | valid samplesheet, paired FASTQ | Controls, sampling plan, independent replicates |
| G1 | fastp/FastQC/MultiQC and PhiX | Residual contamination and library biases; quantify retained pairs |
| G2 | assembly N50, mapped reads, ALE, coverage | Assembly biases, chimeras, alternatives where justified |
| G3 | gene quality, taxonomy conflict matrix | Prodigal unsuitable for euk multi-exon proteins; viral frames |
| G4 | bin CheckM2, GUNC, GTDB, rRNA/tRNA | MIMAG quality, lineage-appropriate QC; manual review |
| G5 | score/resolution for virus+SSU+euk | CheckV quality, PR2 18S tree; avoid species overcalls |
| G6 | KO/HMM orthology and context | Metabolic pathway completeness/directionality, key marker trees |
| G7 | data-derived abundances | Read-vs-contig-vs-gene-vs-MAG denominators, multiple mapping |
| G8 | statistical models (when replicated) | FDR, batch confounding, compositionality, independent units |
| G9 | DNA claims labelled as potential | RNA/protein/assays required for expression/flux conclusions |
| G10 | provenance, manifests, checksums, version pins | MIxS/MIMAG/MIUViG/FAIR and deposition approvals |

Acceptance thresholds should be literature- and dataset-specific; avoid unvalidated universal N50, identity, e-value or MAG completeness cutoffs. Annotate accepted thresholds and rationale in methods before running.
