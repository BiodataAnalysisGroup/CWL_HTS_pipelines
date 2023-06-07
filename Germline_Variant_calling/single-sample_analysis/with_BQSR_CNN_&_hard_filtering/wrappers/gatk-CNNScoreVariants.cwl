cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, CNNScoreVariants]

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
      prefix: -R
      position: 1
      shellQuote: false
  variant:
    type: File
    # secondaryFiles:
    #     - .tbi
    inputBinding:
      prefix: -V
      position: 2
      shellQuote: false
  aligned_reads:
    type: File?
    secondaryFiles:
      - .bai
    inputBinding: 
      prefix: -I
      position: 3
      shellQuote: false
  tensor_type: 
    type: string?
    inputBinding:
      prefix: -tensor-type
      position: 4
      shellQuote: false
  architecture:
    type: File?
    inputBinding:
      prefix: -architecture
      position: 5
  weights:
    type: File?
    inputBinding:
      prefix: -weights
      position: 6
  OUTPUT:
    type: string
    inputBinding:
      prefix: -O
      position: 100
      shellQuote: false
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.OUTPUT)