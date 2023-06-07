cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

hints:
- class: DockerRequirement
  dockerPull: biodataanalysisgroup/deseq2-rscript:v1.0

baseCommand: ["Rscript", "--vanilla"]

inputs:
    rscript: 
        type: string
        default: /home/analysis/DESeq2.R
        inputBinding:
            position: 1
    count_matrix:
        type: File
        inputBinding:
            position: 2
            prefix: --count-matrix
        doc: TSV/CSV file containing RNA-Seq read counts matrix.
    metadata:
        type: File
        inputBinding:
            position: 3
            prefix: --metadata
        doc: TSV/CSV file with metadata on RNA-Seq samples, including the phenotype of interest.
    samples:
        type: string
        inputBinding: 
            position: 4
            prefix: --samples
        doc: Name of the column with sample names/IDs in TSV/CSV file with metadata.
    design:
        type: string
        inputBinding:
            position: 5
            prefix: --design
        doc: Formula-based design for the differential expression analysis.
    min_sum_of_reads:
        type: int?
        default: 10
        inputBinding:
            position: 6
            prefix: --min-sum-of-reads
        doc: Initial minimum threshold for filtering out genes with very low counts (based on their total sum of reads).
    reference_level:
        type: string?
        inputBinding:
            position: 7
            prefix: --reference-level
        doc: "Passes string values, not integers, to the R script to determine the reference (level) label for the phenotype of interest."
    phenotype: 
        type: string?
        inputBinding: 
            position: 8
            prefix: --phenotype
        doc: |
         Variable (phenotype) of interest based on which to perform differential expression analysis. 
         If not defined then the last column of the table passed to the --metadata argument will be used.
    contrast:
        type: boolean?
        default: false
        inputBinding:
            position: 9
            prefix: --contrast
        doc: Specify numerator and denomimator for constrast during DE analysis.
    numerator:
        type: string?
        inputBinding:
            position: 10
            prefix: --numerator
        doc: Phenotype label to be used as the numerator during DE analysis. Requires --contrast to be set to true.
    denominator:
        type: string?
        inputBinding:
            position: 11
            prefix: --denominator
        doc: Phenotype label to be used as the denomimator during DE analysis. Requires --contrast to be set to true.
    lfcThreshold: 
        type: float
        default: 0
        inputBinding: 
            position: 12
            prefix: --lfcThreshold
        doc: Log Fold Change (LFC) threshold for genes to be considered differentially expressed.
    pAdjustMethod:
        type: string
        default: BH
        inputBinding:
            position: 13
            prefix: --pAdjustMethod
        doc: Multiple testing correction method.
    alpha: 
        type: float
        default: 0.1
        inputBinding: 
            position: 14
            prefix: --alpha
        doc: Qvalue threshold for genes to be considered differentially expressed.
    parallelization:
        type: boolean?
        default: false
        inputBinding:
            position: 15
            prefix: --parallelization
        doc: |
         If set to true, parallel execution of DESeq2::DESeq() and DESeq2::results() using BiocParallel. 
         The number of cores for the register function are defined by the --cores argument.
    cores:
        type: int?
        inputBinding:
            position: 16
            prefix: --cores
        doc: Number of cores to run in parallel. Requires --parallelization to be set to true.
    transformation: 
        type: string?
        default: vst
        inputBinding:
            position: 17
            prefix: --transformation
        doc: |
         Transformation method for read counts for other downstream analyses – e.g. for visualization or clustering. 
         Available options are vst (Variance stabilizing transformation) and rlog (Regularized log transformation).
    blind:
        type: boolean?
        default: false
        inputBinding:
            position: 18
            prefix: --blind
        doc: |
         It uses the design formula to calculate the within-group variability (if blind=FALSE) or the across-all-samples variability (if blind=TRUE).
    hypothesis:
        type: string
        default: Wald
        inputBinding: 
            position: 19
            prefix: --hypothesis
        doc: "Type of hypothesis testing: Wald or Likelihood Ratio Test (LRT)"
    reduced:
        type: string?
        inputBinding:
            position: 20
            prefix: --reduced
        doc: Formula string with variable(s) to remove from and compare with the initial DE design during LRT testing. It is required if --hypothesis is set to LRT.
    hidden_batch_effects:
        type: boolean?
        default: false
        inputBinding:
            position: 21
            prefix: --hidden-batch-effects
        doc: Automatically perform analysis for hidden batch effects using SVA or RUV methods for a pre-determined number of factors.
    hidden_batch_row_means:
        type: int?
        default: 10
        inputBinding:
            position: 22
            prefix: --hidden-batch-rowMeans
        doc: The average minimum threshold for each gene to be considered for hidden batch effects analysis.
    hidden_batch_method:
        type: string?
        default: SVA
        inputBinding:
            position: 23
            prefix: --hidden-batch-method
        doc: Use the SVA or RUV method for hidden batch effects analysis.
    variables:
        type: int?
        inputBinding: 
            position: 24
            prefix: --hidden-batch-variables
        doc: Number of factors for which hidden batch effects analysis will be performed.
    
outputs:
    deseq2_de_results:
        type: File
        outputBinding:
            glob: "DESeq2_results.tsv"
    deseq2_dds_object:
        type: File
        outputBinding:
            glob: "DESeq2_DESeqDataSet.rda"
    deseq2_res_lfcShrink_object:
        type: File
        outputBinding:
            glob: "DESeq2_lfcShrink.rda"
    deseq2_transformed_object:
        type: File?
        outputBinding:
            glob: "DESeq2_DESeqTransform.rda"
