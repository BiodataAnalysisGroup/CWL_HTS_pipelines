cwlVersion: v1.0
class: CommandLineTool

requirements:
    - class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastx_toolkit:v0.0.14

arguments:
  - position: 1
    shellQuote: false
    valueFrom: $(inputs.command)
  - prefix: -o
    position: 6
    valueFrom: $(inputs.input_file.basename.split(".fq")[0] + "_trim.fq")

baseCommand: ["/usr/local/bin/fastx_trimmer"]
inputs:
  command:
    type: string?
    default: "/usr/local/bin/fastx_trimmer"
  first_base_to_keep:
    type: int?
    inputBinding:
      prefix: -f
      position: 2
    doc: |
      Default is first base.
  last_base_to_keep:
    type: int?
    inputBinding:
      prefix: -l
      position: 3
    doc: |
      Default is the entire read.
  input_file:
    type: File
    inputBinding:
      prefix: -i
      position: 4
    doc: |
      Input fastq file.
  compressed_output:
    type: boolean?
    default: false
    inputBinding:
      prefix: -z
      position: 5

outputs:
  output:
    type: File?
    outputBinding:
      glob: $(inputs.input_file.basename.split(".fq")[0] + "_trim.fq")
