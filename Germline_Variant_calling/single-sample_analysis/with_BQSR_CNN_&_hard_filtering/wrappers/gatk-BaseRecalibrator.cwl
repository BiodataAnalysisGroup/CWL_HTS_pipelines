cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, BaseRecalibrator]

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
  known_sites_1:
    type: File
    secondaryFiles:
      - .tbi
    inputBinding:
      position: 3
      prefix: --known-sites
  known_sites_2:
    type: File?
    secondaryFiles:
      - .tbi
    inputBinding:
      position: 3
      prefix: --known-sites
  known_sites_3:
    type: File?
    secondaryFiles:
      - .tbi
    inputBinding:
      position: 3
      prefix: --known-sites
  known_sites_4:
    type: File?
    secondaryFiles:
      - .tbi
    inputBinding:
      position: 3
      prefix: --known-sites
  known_sites_5:
    type: File?
    secondaryFiles:
      - .tbi
    inputBinding:
      position: 3
      prefix: --known-sites
  intervals: 
    type: File?
    inputBinding: 
      prefix: -L
      position: 4
    label: One or more genomic intervals over which to operate
  interval_padding:
    type: int?
    inputBinding:
      prefix: --interval-padding
      position: 5
    label: Amount of padding (in bp) to add to each interval you are including.
  OUTPUT:
    type: string

arguments:
  - prefix: -O
    valueFrom: $(inputs.OUTPUT + "." + inputs.intervals.basename.split("-")[0] + ".table")
    position: 100

outputs:
  output_table: 
    type: File
    outputBinding:
      glob: $(inputs.OUTPUT + "." + inputs.intervals.basename.split("-")[0] + ".table") 
