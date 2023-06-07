cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: biodataanalysisgroup/chipqc-diffbind-rscripts:v1.0

baseCommand: ["Rscript", "--vanilla"]

inputs:
    rscript: 
        type: string
        default: /home/analysis/DiffBind.R
        inputBinding:
            position: 1
    treatmentBAM: 
        type: File[]
        secondaryFiles: 
            - .bai
        inputBinding: 
            position: 2
            prefix: --treatment-bam
    controlBAM: 
        type: File[]?
        secondaryFiles: 
            - .bai
        inputBinding: 
            position: 3
            prefix: --control-bam
    peaks: 
        type: File[]
        inputBinding:
            position: 4
            prefix: --peak-files
    peakCaller: 
        type: string
        default: bed
        inputBinding: 
            position: 5
            prefix: --peak-caller
    metadata:
        type: File
        inputBinding:
            position: 6
            prefix: --metadata
    consensus:
        type: string[]?
        inputBinding: 
            position: 7
            prefix: --consensus
        doc: |
         Values from "DiffBind – DBA global constant variables" matching with variables in the metadata table (e.g., "DBA_TISSUE" or "DBA_CONDITION")
    minOverlap:
        type:
        - int
        - float
        default: 2
        inputBinding: 
            position: 8
            prefix: --minOverlap
        doc: |
         Only include peaks in at least this many peaksets in the main binding matrix if basing DBA object on an existing one. 
         If minOverlap is between zero and one, peak will be included from at least this proportion of peaksets.
    summits:
        type: int?
        default: 200
        inputBinding:
            position: 9
            prefix: --summits
    blacklist: 
        type: 
        - string
        - boolean
        - File
        default: true
        inputBinding:
            position: 10
            prefix: --blacklist
        doc: |
         Configure blacklist value (blacklisted genomic regions) as:
         1. auto-detect (set to true to use only --blacklist prefix),
         2. DiffBind - DBA constant variable (e.g., "DBA_BLACKLIST_HG38"), or
         3. Use BED file containing custom regions (e.g., hg38-blacklist.v2.bed ENCODE file)
    greylist: 
        type: 
        - string
        - boolean
        - File
        default: false
        inputBinding:
            position: 11
            prefix: --greylist
        doc: |
         Configure greylist value as:
         1. auto-detect (set to true to use only --greylist prefix),
         2. DiffBind - DBA constant variable (e.g., "DBA_BLACKLIST_HG38") which refers to a BSgenome object, or
         3. Use BED file containing custom regions (will be converted to GRanges object)
    cores: 
        type: int?
        inputBinding: 
            position: 12
            prefix: --cores
    fragmentSize:
        type: int?
        default: 125
        inputBinding:
            position: 13
            prefix: --fragmentSize
    bParallel:
        type: boolean?
        default: true
        inputBinding: 
            position: 14
            prefix: --bParallel
    normalization: 
        type: string
        default: "DBA_NORM_DEFAULT"
        inputBinding: 
            position: 15
            prefix: --normalization
        doc: |
         Use only DiffBind – DBA global constant variables
    library:
        type: string
        default: "DBA_LIBSIZE_DEFAULT"
        inputBinding: 
            position: 16
            prefix: --library
        doc: |
         Use only DiffBind – DBA global constant variables
    background: 
        type: boolean?
        inputBinding: 
            position: 17
            prefix: --background
        doc: |
         Apply or not background normalization
    spikein: 
        type: 
        - boolean?
        - File?
        default: false
        inputBinding: 
            position: 18
            prefix: --spikein
        doc: |
         Assign spikein using tracks (one for each sample in metadata table) (value: true) or for parallel factor normalization (File).
         Requires further development, operates with the default 'false' for this version.
    design: 
        type: string
        default: "~Condition"
        inputBinding: 
            position: 19
            prefix: --design
        doc: |
         Use a valid formula design for differential binding comparison between samples. 
         Only columns from the metadata table should be used (e.g., "Condition" or "Tissue"). 
         Multifactor designs can be provided (e.g., "~Tissue + Condition").
    contrast:
        type: string[]?
        inputBinding: 
            position: 20
            prefix: --contrast
        doc: |
         Character vector of length three, to define specific contrasts.
    reorderMeta_factor:
        type: string[]?
        inputBinding: 
            position: 21
            prefix: --reorderMeta-factor
        doc: |
         Metadata table column name to connect 'reorderMeta_value' and assign reference values (denominator in dba.analysis).
    reorderMeta_value:
        type: string[]?
        inputBinding: 
            position: 22
            prefix: --reorderMeta-value
        doc: |
         Value from the contents of the respective metadata table column ('reorderMeta_factor') that will be assigned as reference (denominator in dba.analysis).
    retrieve_consensus:
        type: boolean?
        default: true
        inputBinding: 
            position: 23
            prefix: --retrieve-consensus
    diffbind_results_filename:
        type: string
        inputBinding: 
            position: 24
            prefix: --diffbind-results
    correlation_heatmap_filename: 
        type: string
        inputBinding:
            position: 25
            prefix: --correlation-heatmap
    diffbind_consensus_filename:
        type: string
        inputBinding:
            position: 26
            prefix: --diffbind-consensus
    low_read_count_filter:
        type: int
        default: 1
        inputBinding:
            position: 27
            prefix: --low-read-count-filter
    diffbind_filterFun:
        type: string
        default: "max"
        inputBinding:
            position: 28
            prefix: --diffbind-filterFun

outputs:
    diffbind_results:
        type: File
        outputBinding:
            glob: $(inputs.diffbind_results_filename) 
    correlation_heatmap:
        type: File?
        outputBinding:
            glob: $(inputs.correlation_heatmap_filename) 
    diffbind_consensus:
        type: File?
        outputBinding:
            glob: $(inputs.diffbind_consensus_filename) 
    diffbind_normalized_counts:
        type: File?
        outputBinding:
            glob: $("DiffBind_normalized_counts.csv")
    diffbind_dba_object:
        type: File?
        outputBinding:
            glob: $("DiffBind_DBA_object.rda")