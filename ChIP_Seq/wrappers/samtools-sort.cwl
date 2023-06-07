cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/samtools:1.14--hb421002_0"

baseCommand: [samtools, sort]

inputs:
  compression_level:
    type: int?
    inputBinding:
      prefix: -l
    doc: |
      Set compression level, from 0 (uncompressed) to 9 (best)
  threads:
    type: int?
    inputBinding:
      prefix: -@

    doc: Set number of sorting and compression threads [1]
  memory:
    type: string?
    inputBinding:
      prefix: -m
    doc: |
      Set maximum memory per thread; suffix K/M/G recognized [768M]
  input:
    type: File
    inputBinding:
      position: 1
    doc: Input bam file.
  output_name:
    type: string?
    inputBinding:
      position: 2
      prefix: -o
    doc: Desired output filename.
  sort_by_name:
    type: boolean?
    default: false
    inputBinding:
      prefix: -n
    doc: Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
outputs:
  sorted:
    type: File
    outputBinding:
      glob: $(inputs.output_name)