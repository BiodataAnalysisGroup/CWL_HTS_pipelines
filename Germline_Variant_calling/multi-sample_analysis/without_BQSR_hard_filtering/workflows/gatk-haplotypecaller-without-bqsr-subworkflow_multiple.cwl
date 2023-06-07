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
    bqsr_intervals:
        type: File[]
    bqsr_hc_native_pairHMM_threads:
        type: int?

outputs:
    o_gatk_HaplotypeCaller:
        type: File[]
        outputSource: gatk_HaplotypeCaller/output
    o_gatk_MergeVCFs:
        type: File
        outputSource: gatk_MergeVCFs/output

steps:
    ### GATK - HaplotypeCaller ###
    gatk_HaplotypeCaller:
        run: ../wrappers/gatk-HaplotypeCaller_multiple_v1.cwl 
        # scatterMethod: dotproduct
        scatter:
        - intervals
        in:
            INPUT: bqsr_INPUT
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
