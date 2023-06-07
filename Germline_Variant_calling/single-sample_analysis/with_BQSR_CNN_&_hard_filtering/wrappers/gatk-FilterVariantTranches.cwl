cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0 

baseCommand: [gatk, FilterVariantTranches]

inputs:
  variant:
    type: File
    # secondaryFiles:
    #     - .tbi
    inputBinding:
      prefix: -V
      position: 1
      shellQuote: false
  resource_1:
    type: File
    secondaryFiles: 
        - .tbi
    inputBinding:
      prefix: --resource
      position: 2
      shellQuote: false
  resource_2:
    type: File?
    secondaryFiles: 
        - .tbi
    inputBinding:
      prefix: --resource
      position: 2
      shellQuote: false    
  resource_3:
    type: File?
    secondaryFiles: 
        - .tbi
    inputBinding:
      prefix: --resource
      position: 2
      shellQuote: false
  info_key:
    type: string
    default: CNN_2D
    inputBinding:
      prefix: --info-key
      position: 3
      shellQuote: false
  SNP:
    type: float
    default: 99.95
    inputBinding: 
      prefix: --snp-tranche
      position: 4
      shellQuote: false
  INDEL:
    type: float
    default: 99.4
    inputBinding: 
      prefix: --indel-tranche
      position: 5
      shellQuote: false
  invalidate_previous_filters:
    type: boolean?
    inputBinding:
      prefix: --invalidate-previous-filters
      position: 6
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
    outputBinding:
      glob: $(inputs.OUTPUT)