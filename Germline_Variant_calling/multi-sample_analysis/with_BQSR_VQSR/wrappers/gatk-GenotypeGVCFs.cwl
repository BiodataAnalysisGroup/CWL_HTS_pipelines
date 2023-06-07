cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, GenotypeGVCFs]

inputs:
    reference:
        type: File
        secondaryFiles:
        - .amb
        - .ann
        - .bwt
        - .pac
        - .sa
        - .fai
        - ^.dict
        inputBinding: 
            position: 1
            prefix: -R
    gvcf_input:
        type: File
        secondaryFiles:
        - .tbi
        inputBinding:
            position: 2
            prefix: -V
    vcf_output:
        type: string
        default: "cohort.genotyped.vcf.gz"
        inputBinding:
            position: 3
            prefix: -O

outputs:
    output:
        type: File
        outputBinding:
            glob: $(inputs.vcf_output)
        secondaryFiles: 
        - .tbi