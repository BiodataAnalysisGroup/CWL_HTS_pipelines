cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, SelectVariants]

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
    secondaryFiles: 
        - .tbi
    inputBinding:
      prefix: -V
      position: 2
      shellQuote: false
  exclude_filter:
    type: boolean?
    inputBinding:
      prefix: --exclude-filtered
      position: 3
      shellQuote: false
  variant_type:
    type: string?
    inputBinding:
      prefix: --select-type-to-include
      position: 4
      shellQuote: false
  rm_unused_alt:
    type: boolean?
    inputBinding:
      prefix: --remove-unused-alternates
      position: 5
      shellQuote: false
  OUTPUT:
    type: string
    inputBinding:
      prefix: -O
      position: 100
      shellQuote: false

outputs:
  output:
    type: File
    secondaryFiles: 
        - .tbi
    outputBinding:
      glob: $(inputs.OUTPUT)