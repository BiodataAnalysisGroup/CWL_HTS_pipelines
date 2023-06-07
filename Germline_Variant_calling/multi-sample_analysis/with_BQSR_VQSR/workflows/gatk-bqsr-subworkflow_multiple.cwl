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
    bqsr_hc_native_pairHMM_threads:
        type: int?
    bqsr_hc_java_options:
        type: string?

outputs:
    o_gatk_baserecalibrator_table:
        type: File[]
        outputSource: gatk_baserecalibrator/output_table
    o_gatk_applybqsr:
        type: File[]
        outputSource: gatk_applybqsr/output
    o_gatk_HaplotypeCaller:
        type: File[]
        outputSource: gatk_HaplotypeCaller/output
    o_gatk_MergeVCFs:
        type: File
        outputSource: gatk_MergeVCFs/output

steps:
    ### GATK - BaseRecalibrator ###
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
    ### GATK - ApplyBQSR ###
    gatk_applybqsr:
        run: ../wrappers/gatk-ApplyBQSR_v2.cwl
        scatterMethod: dotproduct
        scatter: 
        - intervals
        - bqsr_recal_file
        in:
            INPUT: bqsr_INPUT
            reference: bqsr_reference
            intervals: bqsr_intervals
            bqsr_recal_file: gatk_baserecalibrator/output_table
            OUTPUT: bqsr_OUTPUT
        out: [output]
    ### GATK - HaplotypeCaller ###
    gatk_HaplotypeCaller:
        run: ../wrappers/gatk-HaplotypeCaller_multiple_v1.cwl 
        scatterMethod: dotproduct
        scatter:
        - INPUT
        - intervals
        in:
            java_options: bqsr_hc_java_options
            INPUT: gatk_applybqsr/output
            OUTPUT: bqsr_OUTPUT
            reference: bqsr_reference
            intervals: bqsr_intervals
            native_pairHMM_threads: bqsr_hc_native_pairHMM_threads
        out: [output]
    ### GATK - MergeVCFs ###
    gatk_MergeVCFs:
        run: ../wrappers/gatk-MergeVCFs.cwl 
        in: 
            INPUT: gatk_HaplotypeCaller/output
            OUTPUT: bqsr_OUTPUT
        out: [output]
