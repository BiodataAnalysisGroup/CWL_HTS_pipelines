cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, ApplyBQSR]

inputs:
  INPUT: 
    type: File
    inputBinding:
      position: 1
      prefix: -I
      shellQuote: false
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
      prefix: -R
      position: 2
      shellQuote: false
  bqsr_recal_file:
    type: File
    inputBinding:
      prefix: --bqsr-recal-file
      position: 3
      shellQuote: false
  OUTPUT:
    type: string
    inputBinding:
      prefix: -O
      position: 4
      shellQuote: false

outputs:
  output: 
    type: File
    outputBinding:
      glob: $(inputs.OUTPUT)
