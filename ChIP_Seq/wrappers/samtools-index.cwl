cwlVersion: v1.0
class: CommandLineTool

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.alignments) ]
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/samtools:1.14--hb421002_0"

baseCommand: [samtools, index]

inputs:
  alignments:
    type: File
    inputBinding:
      valueFrom: $(self.basename)
    label: Input bam file.

outputs:
  alignments_with_index:
    type: File
    secondaryFiles: .bai
    outputBinding:
      glob: $(inputs.alignments.basename)
