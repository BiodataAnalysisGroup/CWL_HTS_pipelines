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
    secondaryFiles:
      - .bai
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
  OUTPUT:
    type: string

arguments:
- prefix: -O
  valueFrom: $(inputs.OUTPUT + "." + inputs.intervals.basename.split("-")[0] + ".bqsr.bam")

outputs:
  output: 
    type: File
    secondaryFiles:
      - ^.bai
    outputBinding:
      glob: $(inputs.OUTPUT + "." + inputs.intervals.basename.split("-")[0] + ".bqsr.bam")
