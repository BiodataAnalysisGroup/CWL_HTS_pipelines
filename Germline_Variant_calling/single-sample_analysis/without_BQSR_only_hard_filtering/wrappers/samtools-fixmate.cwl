cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.14--hb421002_0"

baseCommand: [samtools, fixmate]

inputs:
  threads:
    type: int
    default: 16
    inputBinding:
      prefix: -@
      position: 1
  output_format:
    type: string
    default: bam
    inputBinding:
      prefix: -O
      position: 1
  r: 
    type: float?
    inputBinding:
      prefix: -r
      position: 2
    doc: Remove secondary and unmapped reads.
  p: 
    type: float?
    inputBinding:
      prefix: -p
      position: 3
    doc: Disable FR proper pair check.
  c: 
    type: float?
    inputBinding:
      prefix: -c
      position: 4
    doc: Add template cigar ct tag.
  m: 
    type: float?
    inputBinding:
      prefix: -m
      position: 5
    doc: Add ms (mate score) tags. These are used by markdup to select the best reads to keep.
  uncompressed_output:
    type: float?
    inputBinding:
      prefix: -u
      position: 6
    doc: Output uncompressed BAM or CRAM.
  no_PG:
    type: float?
    inputBinding:
      prefix: --no-PG
      position: 7
    doc: Do not add a @PG line to the header of the output file.
  input_file:
    type: File
    inputBinding:
      position: 100
  output_file_name:
    type: string
    inputBinding:
      position: 101
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)