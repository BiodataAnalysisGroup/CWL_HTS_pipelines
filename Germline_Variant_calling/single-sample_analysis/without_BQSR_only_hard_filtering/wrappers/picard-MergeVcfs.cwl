cwlVersion: v1.0
class: CommandLineTool

requirements:
    DockerRequirement:
        dockerPull: "quay.io/biocontainers/picard:2.26.7--hdfd78af_0"

baseCommand: [picard, MergeVcfs]

inputs:
  INPUT_SNP:
    type: File
    inputBinding:
      separate: false
      prefix: I=
      position: 1
  INPUT_INDEL:
    type: File
    inputBinding:
        separate: false
        prefix: I=
        position: 2
  OUTPUT:
    type: string
    inputBinding:
      separate: false
      prefix: O=
      position: 3

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.OUTPUT)
