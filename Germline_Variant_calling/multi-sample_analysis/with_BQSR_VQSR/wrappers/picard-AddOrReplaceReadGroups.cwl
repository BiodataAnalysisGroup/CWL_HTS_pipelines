cwlVersion: v1.0
class: CommandLineTool

requirements:
    DockerRequirement:
        dockerPull: "quay.io/biocontainers/picard:2.26.7--hdfd78af_0"

baseCommand: [picard, AddOrReplaceReadGroups]

inputs:
  INPUT:
    type: File
    inputBinding:
      separate: false
      prefix: I=
      position: 3
  OUTPUT:
    type: string?
    inputBinding:
      separate: false
      prefix: O=
      position: 4
  rgid:
    type: string?
    inputBinding:
      separate: false
      prefix: RGID=
      position: 5
    doc: Read-Group ID
  rglb:
    type: string?
    inputBinding:
      separate: false
      prefix: RGLB=
      position: 5
  rgpl:
    type: string?
    inputBinding:
      separate: false
      prefix: RGPL=
      position: 5
  rgpu:
    type: string?
    inputBinding:
      separate: false
      prefix: RGPU=
      position: 5
  rgsm:
    type: string?
    inputBinding:
      separate: false
      prefix: RGSM=
      position: 5

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.OUTPUT)
