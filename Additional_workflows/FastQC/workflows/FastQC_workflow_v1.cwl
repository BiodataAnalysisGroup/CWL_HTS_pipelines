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
    # Options
    fastqc_command:
        type: string?
    cp_command: 
        type: string?
    
outputs:
    ### FASTQC outputs ###
    # HTML
    o_fastqc_html:
        type: File[]?
        outputSource: rename_fastqc_html/renamed_file
    # ZIP
    o_fastqc_zip:
        type: Directory?
        outputSource: cp_fastqc_zip/output_dir

steps:
    ### Retrieve FASTQ files from target directory (ExpressionTool) ###
    get_raw_files:
        run: ../wrappers/get-raw-files.cwl
        in:
            DIRECTORY: raw_files_directory
        out: [raw_files]
    ### FASTQC - FASTQ file quality control ###
    fastqc:
        run: ../wrappers/fastqc.cwl
        scatter: 
            - input_files
        in:
            command: fastqc_command
            input_files: get_raw_files/raw_files
        out: 
            - html_file
            - zipped_file 
    ### Copy FASTQC output files - ZIP ###
    cp_fastqc_zip:
        run: ../wrappers/cp.cwl
        in:
            command: cp_command
            input_files: fastqc/zipped_file
            outputdir:
                valueFrom: $( "fastqc_zip" )
        out:
            - output_dir
    ### Rename FASTQC output files - HTML ###
    rename_fastqc_html:
        run: ../wrappers/rename.cwl
        in:
            input_file: fastqc/html_file
            run_type:
                valueFrom: $( "fastqc_" )
        out:
            - renamed_file
        scatter:
            - input_file
