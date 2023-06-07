cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "biocontainers/bcftools:v1.5_cv3"

baseCommand: [bcftools, view]
inputs:
  input: 
    type: File
    inputBinding: 
        position: 100
  threads: 
    type: int?
    inputBinding: 
        prefix: --threads
        position: 1
  include: 
    type: string?
    inputBinding: 
        prefix: --include
        position: 2
  exclude: 
    type: string?
    inputBinding: 
        prefix: --exclude
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
        prefix: --output-file
  
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_name)