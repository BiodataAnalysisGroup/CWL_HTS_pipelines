cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.14--hb421002_0"

baseCommand: [samtools, fixmate]

arguments: ["-m"] # Add ms (mate score) tags. These are used by markdup to select the best reads to keep.

inputs:
  threads:
    type: int
    default: 8
    inputBinding:
      prefix: -@
      position: 1
  output_format:
    type: string
    default: bam
    inputBinding:
      prefix: -O
      position: 1
  input_file:
    type: File
    inputBinding:
      position: 2
  output_file_name:
    type: string
    inputBinding:
      position: 3
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)