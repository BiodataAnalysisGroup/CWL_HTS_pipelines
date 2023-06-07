class: Workflow

cwlVersion: v1.0

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement

inputs:
    bqsr_INPUT:
        type: File
        secondaryFiles: 
            - .bai
    bqsr_OUTPUT:
        type: string
    bqsr_reference:
        type: File
        secondaryFiles:
            - .amb
            - .ann
            - .bwt
            - .pac
            - .sa
            - .fai
            - ^.dict
    bqsr_known_sites_1:
        type: File
        secondaryFiles: 
            - .tbi
    bqsr_known_sites_2:
        type: File?
        secondaryFiles: 
            - .tbi
    bqsr_known_sites_3:
        type: File?
        secondaryFiles: 
            - .tbi
    bqsr_intervals:
        type: File[]
    bqsr_interval_padding:
        type: int?

outputs:
    o_gatk_baserecalibrator_table:
        type: File[]
        outputSource: gatk_baserecalibrator/output_table
    o_gatk_gatherbqsrreports:
        type: File
        outputSource: gatk_gatherbqsrreports/output

steps:
    gatk_baserecalibrator:
        run: ../wrappers/gatk-BaseRecalibrator.cwl
        scatter: 
        - intervals
        in: 
            INPUT: bqsr_INPUT
            OUTPUT: bqsr_OUTPUT
            reference: bqsr_reference
            known_sites_1: bqsr_known_sites_1
            known_sites_2: bqsr_known_sites_2
            known_sites_3: bqsr_known_sites_3
            intervals: bqsr_intervals
            interval_padding: bqsr_interval_padding
        out: [output_table]
    # gatk GatherBQSRReports
    gatk_gatherbqsrreports:
        run: ../wrappers/gatk-GatherBQSRReports.cwl
        in: 
            bqsr_tables: gatk_baserecalibrator/output_table
            output_file: bqsr_OUTPUT
        out: [output]
