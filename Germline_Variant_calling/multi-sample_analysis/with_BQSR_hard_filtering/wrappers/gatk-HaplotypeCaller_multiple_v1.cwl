cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, HaplotypeCaller]

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
  INPUT:
    type: File
    inputBinding:
      prefix: -I
      position: 2
      shellQuote: false
    secondaryFiles:
      - ^.bai
  native_pairHMM_threads:
    type: int?
    inputBinding:
      prefix: --native-pair-hmm-threads
      position: 3
  ploidy:
    type: int?
    inputBinding:
      prefix: --sample-ploidy
      position: 4
  intervals:
    type: File?
    inputBinding:
      prefix: -L
      position: 5
      shellQuote: false
  exclude_intervals:
    type: File?
    inputBinding:
      prefix: -XL
      position: 6
      shellQuote: false
  max_mnp_distance:
    type: int
    default: 0
    inputBinding:
      position: 7
      prefix: --max-mnp-distance
      shellQuote: false
  erc: 
    type: string
    default: GVCF
    inputBinding:
        prefix: -ERC
        position: 8
        shellQuote: false
  OUTPUT:
    type: string
    # inputBinding:
    #   prefix: -O
    #   position: 100
    #   shellQuote: false

arguments:
- prefix: -O
  valueFrom: $(inputs.OUTPUT + "." + inputs.intervals.basename.split("-")[0] + ".raw.g.vcf.gz")

outputs:
  output:
    type: File
    secondaryFiles: 
    - .tbi
    outputBinding:
      glob: $(inputs.OUTPUT + "." + inputs.intervals.basename.split("-")[0] + ".raw.g.vcf.gz")