cwlVersion: v1.0
class: Workflow

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
    ### Trimmomatic ###
    # SE
    trimmomatic_se_threads: 
        type: int?
    trimmomatic_se_illuminaClip:
        type: string?
    trimmomatic_se_slidingWindow:
        type: string?
    trimmomatic_se_leading:
        type: int?
    trimmomatic_se_trailing:
        type: int?
    trimmomatic_se_minlen:
        type: int?
    # PE
    trimmomatic_pe_threads: 
        type: int?
    trimmomatic_pe_illuminaClip:
        type: string?
    trimmomatic_pe_slidingWindow:
        type: string?
    trimmomatic_pe_leading:
        type: int?
    trimmomatic_pe_trailing:
        type: int?
    trimmomatic_pe_minlen:
        type: int?

outputs:
    # Trimmomatic outputs
    o_trimmomatic_single_end_stderr:
        type: File[]
        outputSource: trimmomatic_single_end/stderr_log
    o_trimmomatic_single_end_fastq:
        type: File[]
        outputSource: trimmomatic_single_end/outFastq
    o_trimmomatic_paired_end_stderr:
        type: File[]
        outputSource: trimmomatic_paired_end/stderr_log
    o_trimmomatic_paired_end_fwd_paired:
        type: File[]
        outputSource: trimmomatic_paired_end/outFastq_fwd_paired
    o_trimmomatic_paired_end_fwd_unpaired:
        type: File[]
        outputSource: trimmomatic_paired_end/outFastq_fwd_unpaired
    o_trimmomatic_paired_end_rev_paired:
        type: File[]
        outputSource: trimmomatic_paired_end/outFastq_rev_paired
    o_trimmomatic_paired_end_rev_unpaired:
        type: File[]
        outputSource: trimmomatic_paired_end/outFastq_rev_unpaired

steps:
    # FASTQC QUALITY CHECK - TRIMMMING - MAPPING - FILTER MAPPED READS
    ### Separate files - Generate file names (ExpressionTool) ###
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
            - trimmomatic_command_single
            - trimmomatic_command_paired
    ### Trimmomatic ###
    trimmomatic_single_end:
        run: ../wrappers/Trimmomatic_SE.cwl
        scatter:
            - input_fastq
        in: 
            input_fastq: split_single_paired/single_files
            trimm_se_threads: trimmomatic_se_threads
            command: split_single_paired/trimmomatic_command_single
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            illuminaClip: trimmomatic_se_illuminaClip
            slidingWindow: trimmomatic_se_slidingWindow
            leading: trimmomatic_se_leading
            trailing: trimmomatic_se_trailing
            minlen: trimmomatic_se_minlen
        out:
            - stderr_log
            - outFastq
    trimmomatic_paired_end:
        run: ../wrappers/Trimmomatic_PE.cwl
        scatter:
            - input_fastq_fwd
            - input_fastq_rev
        scatterMethod: dotproduct
        in: 
            input_fastq_fwd: split_single_paired/paired_files_fwd
            input_fastq_rev: split_single_paired/paired_files_rev
            trimm_pe_threads: trimmomatic_pe_threads
            command: split_single_paired/trimmomatic_command_paired
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            file_split_rev: input_file_split_rev
            illuminaClip: trimmomatic_pe_illuminaClip
            slidingWindow: trimmomatic_pe_slidingWindow
            leading: trimmomatic_pe_leading
            trailing: trimmomatic_pe_trailing
            minlen: trimmomatic_pe_minlen
        out:
            - stderr_log
            - outFastq_fwd_paired
            - outFastq_fwd_unpaired
            - outFastq_rev_paired
            - outFastq_rev_unpaired