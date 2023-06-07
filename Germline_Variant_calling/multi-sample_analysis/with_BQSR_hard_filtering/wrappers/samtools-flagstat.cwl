cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/samtools:1.14--hb421002_0"

baseCommand: [samtools, flagstat]

inputs:
  threads:
    type: int?
    default: 16
    inputBinding:
      prefix: -@
      position: 1
  input:
    type: File?
    inputBinding:
      position: 2
    doc: Input bam file.
  output_name: string
stdout: $(inputs.output_name)
outputs:
  output:
    type: stdout