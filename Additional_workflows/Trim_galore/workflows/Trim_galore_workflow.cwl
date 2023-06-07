class: Workflow
cwlVersion: v1.0

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:
    # Directory with fastq(.gz) input files
    raw_files_directory:
        type: Directory
    ### FASTQ file split ###
    input_file_split:
        type: string?
        default: "_R"
    input_file_split_fwd_single:
        type: string?
        default: "R1"
    input_file_split_rev:
        type: string?
        default: "R2"
    ### Trimming options ###
    tg_quality: int
    tg_length: int
    tg_compression: boolean
    tg_do_not_compress: boolean
    tg_trim_suffix: string
    tg_strigency: int

outputs:
    ### Trim_galore outputs ###
    o_trim_galore_single_fq:
        type: File[]
        outputSource: trim_galore_single/trim_galore
    o_trim_galore_single_reports:
        type: File[]
        outputSource: trim_galore_single/trim_galore_report
    o_trim_galore_paired_fq:
        type: File[]
        outputSource: trim_galore_paired/trim_galore
    o_trim_galore_paired_reports:
        type: File[]
        outputSource: trim_galore_paired/trim_galore_report

steps:
    ### Retrieve FASTQ files from target directory (ExpressionTool) ###
    get_raw_files:
        run: ../wrappers/get-raw-files.cwl
        in:
            DIRECTORY: raw_files_directory
        out: [raw_files]
    ### Separate files - Generate file names (ExpressionTool) ###
    split_single_paired:
        run: ../wrappers/split-single-paired_v2.cwl
        in:
            input_files: get_raw_files/raw_files
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            file_split_rev: input_file_split_rev
        out: 
            - single_files
            - paired_files_fwd
            - paired_files_rev
            - trim_galore_command_single
            - trim_galore_command_paired
    ### Trim_galore - FASTQ file trimming ###
    trim_galore_single:
        run: ../wrappers/trim-galore.cwl
        scatter:
            - fq_file_fwd
        in:
            command: split_single_paired/trim_galore_command_single
            fq_file_fwd: split_single_paired/single_files
            length: tg_length
            quality: tg_quality
            compression: tg_compression
            do_not_compress: tg_do_not_compress
            trim_suffix: tg_trim_suffix
            strigency: tg_strigency
            paired: 
                valueFrom: $( false )
        out:
            - trim_galore
            - trim_galore_report
    trim_galore_paired:
        run: ../wrappers/trim-galore.cwl
        scatter:
            - fq_file_fwd
            - fq_file_rev
        scatterMethod: dotproduct
        in:
            command: split_single_paired/trim_galore_command_paired
            fq_file_fwd: split_single_paired/paired_files_fwd
            fq_file_rev: split_single_paired/paired_files_rev
            length: tg_length
            quality: tg_quality
            compression: tg_compression
            do_not_compress: tg_do_not_compress
            trim_suffix: tg_trim_suffix
            strigency: tg_strigency
            paired:
                valueFrom: $( true )
        out: 
            - trim_galore
            - trim_galore_report