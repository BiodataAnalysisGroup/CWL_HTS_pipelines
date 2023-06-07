class: Workflow

cwlVersion: v1.0

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement

inputs:
    hf_INPUT:
        type: File
        secondaryFiles: 
            - .bai
    hf_OUTPUT:
        type: string
    hf_reference:
        type: File
        secondaryFiles:
            - .amb
            - .ann
            - .bwt
            - .pac
            - .sa
            - .fai
            - ^.dict
    hf_intervals:
        type: File[]
    hf_hc_native_pairHMM_threads:
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
            INPUT: hf_INPUT
            OUTPUT: hf_OUTPUT
            reference: hf_reference
            intervals: hf_intervals
            native_pairHMM_threads: hf_hc_native_pairHMM_threads
        out: [output]
    ### GATK - MergeVCFs ###
    gatk_MergeVCFs:
        run: ../wrappers/gatk-MergeVCFs.cwl
        in: 
            INPUT: gatk_HaplotypeCaller/output
            OUTPUT: hf_OUTPUT
        out: [output]
