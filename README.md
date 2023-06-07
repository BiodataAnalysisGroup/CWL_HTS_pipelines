# CWL pipelines for the analysis of Next-Generation Sequencing (NGS) data in GenOptics

Pipelines for the asychronous analysis of raw NGS data, incorporating various software tools by using the [Common Workflow Language](https://www.commonwl.org/) (CWL), a standard for describing computational data-analysis workflows. Specifically, they are written as [CWL workflows](https://www.commonwl.org/v1.2/Workflow.html) and include a number of [CWL Command Line tools](https://www.commonwl.org/v1.0/CommandLineTool.html) (wrappers), which make use of [Docker](https://www.docker.com/) containers. 

The following pipelines are currently available and accept both single- and paired-end NGS raw data (FASTQ format):

- RNA-Seq workflow:
    - Trimming of raw reads with *Trim galore*
    - Quality assessment with *FASTQC*
    - Mapping reads with *HISAT2*
    - SAM/BAM files processing with *SAMtools*
    - Assembly of RNA-Seq alignments into potential transcripts with *StringTie*
    - Read counts with *featureCounts*
- Whole Exome-Seq workflow (under development):
    - Trimming of raw reads with *Trim galore* (to be added)
    - Quality assessment with *FASTQC* (to be added)
    - Mapping reads with *BWA-mem*
    - SAM/BAM files processing with *SAMtools*
    - Variant calling with *GATK HaplotypeCaller*
    - Variant filtration and selection with *GATK VariantFiltration* and *SelectVariants*, respectively
    - Variant annotation with *ANNOVAR*
- ChIP-Seq workflow ([CERTH](https://github.com/GenOptics/pipelines/blob/master/ChIP-seq/ChIP-seq_CERTH.md)) (under development):
    - Trimming of raw reads with *Trim galore*
    - Quality assessment with *FASTQC*
    - Mapping reads with *BWA-backtrack*
    - SAM/BAM files processing with *SAMtools*
    - Mark duplicates with *Picard tools*
    - Call peaks with *MACS2*
    - Output processing with *BEDTools*
- ChIP-Seq workflow ([FORTH](https://github.com/GenOptics/pipelines/blob/master/ChIP-seq/ChIP-seq_FORTH_Steps.md)) (under development)
    - Trimming of raw reads with *Trimmomatic*
    - Quality assessment with *FASTQC*
    - Mapping reads with *HISAT2*
    - SAM/BAM files processing with *SAMtools*
    - Remove duplicates with *SAMtools markdup*
    - Calculate quality metrics with *deeptools*
    - Call peaks with *MACS2* 
    - Output processing with *BEDTools*
    - Rank ordering of super-enhancers with *ROSE*

The CWL workflows and wrappers are compatible and can be executed in the [HYPATIA](https://hypatia.athenarc.gr/) Cloud infrastructure which has been developed to support the computational needs of the [ELIXIR-GR community](https://elixir-greece.org/), but also the broader community of life scientists in Greece and abroad. Detailed instructions on how to use HYPATIA can be found [here](https://hypatia.athenarc.gr/index.php?r=site%2Fhelp) and in this [video capture](https://youtu.be/pupWkmkdGhk) of a webinar on ELIXIR-GR's HYPATIA Cloud Infrastructure.
