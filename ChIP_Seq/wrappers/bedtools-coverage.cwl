cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "staphb/bedtools:2.30.0"
- class: InlineJavascriptRequirement

baseCommand: [bedtools, coverage]

inputs:
  reference_file:
    type: File
    inputBinding:
      position: 1
      prefix: -a

  input_file:
    type: File
    inputBinding:
      position: 2
      prefix: -b

  depth_report:
    type: boolean
    default: false
    inputBinding:
      position: 3
      prefix: -d

  count_of_overlaps:
    type: boolean?
    default: true
    inputBinding:
      position: 4
      prefix: -counts
    doc: Only report the count of overlaps, don't compute fraction, etc.

  min_overlap_a: 
    type: float?
    inputBinding:
      position: 5
      prefix: -f
    doc: Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).

  min_overlap_b: 
    type: float?
    inputBinding:
      position: 6
      prefix: -F
    doc: Minimum overlap required as a fraction of B. Default is 1E-9 (i.e., 1bp).

  fraction_overlaps:
    type: float?
    inputBinding:
      position: 7
      prefix: -r
    doc: Require that the fraction of overlap be reciprocal for A and B.

stdout: $(inputs.input_file.nameroot + "_counts.bed")

outputs:
  output:
    type: stdout