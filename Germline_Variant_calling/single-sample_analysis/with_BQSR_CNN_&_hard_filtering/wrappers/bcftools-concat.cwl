cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "biocontainers/bcftools:v1.5_cv3"

baseCommand: [bcftools, concat]
inputs:
  input1: 
    type: File
    inputBinding: 
        position: 100
    secondaryFiles:
      - .tbi
  input2: 
    type: File
    inputBinding: 
        position: 101
    secondaryFiles:
      - .tbi
  threads: 
    type: int?
    inputBinding: 
        prefix: --threads
        position: 1
  allow_overlaps:
    type: boolean
    default: true
    inputBinding:
        prefix: -a
        position: 2
  remove_duplicates:
    type: string?
    inputBinding:
        prefix: --rm-dups
        position: 3
  output_type: 
    type: string?
    default: v
    inputBinding:
        separate: false
        position: 4
        prefix: -O
  output_name:
    type: string
    inputBinding:
        position: 5
        prefix: -o
  
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_name)