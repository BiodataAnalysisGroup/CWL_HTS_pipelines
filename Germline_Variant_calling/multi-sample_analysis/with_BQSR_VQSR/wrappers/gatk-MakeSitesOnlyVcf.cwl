cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, MakeSitesOnlyVcf]

inputs:
    INPUT:
        type: File
        secondaryFiles: 
        - .tbi
        inputBinding:
            prefix: INPUT=
            separate: false
            position: 1
    OUTPUT:
        type: string
        default: "cohort.genotyped.sitesonly.vcf.gz"
        inputBinding:
            prefix: OUTPUT=
            separate: false
            position: 2

outputs:
    output:
        type: File
        secondaryFiles: 
        - .tbi
        outputBinding: 
            glob: $(inputs.OUTPUT)
