cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

arguments:
- position: 1
  shellQuote: false
  valueFrom: gatk
- position: 2
  shellQuote: false
  prefix: --java-options
  valueFrom: $(inputs.java_options)
- position: 3
  shellQuote: false
  valueFrom: HaplotypeCaller
- position: 4
  prefix: -O
  valueFrom: $(inputs.OUTPUT + "." + inputs.intervals.basename.split("-")[0] + ".raw.g.vcf.gz")

inputs:
  java_options:
    type: string
    default: "'-Xmx4g'"
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
      position: 5
      shellQuote: false
  INPUT:
    type: File
    inputBinding:
      prefix: -I
      position: 6
      shellQuote: false
    secondaryFiles:
      - ^.bai
  native_pairHMM_threads:
    type: int?
    inputBinding:
      prefix: --native-pair-hmm-threads
      position: 7
  ploidy:
    type: int?
    inputBinding:
      prefix: --sample-ploidy
      position: 8
  intervals:
    type: File?
    inputBinding:
      prefix: -L
      position: 9
      shellQuote: false
  exclude_intervals:
    type: File?
    inputBinding:
      prefix: -XL
      position: 10
      shellQuote: false
  max_mnp_distance:
    type: int
    default: 0
    inputBinding:
      position: 11
      prefix: --max-mnp-distance
      shellQuote: false
  erc: 
    type: string
    default: GVCF
    inputBinding:
        prefix: -ERC
        position: 12
        shellQuote: false
  OUTPUT:
    type: string

outputs:
  output:
    type: File
    secondaryFiles: 
    - .tbi
    outputBinding:
      glob: $(inputs.OUTPUT + "." + inputs.intervals.basename.split("-")[0] + ".raw.g.vcf.gz")