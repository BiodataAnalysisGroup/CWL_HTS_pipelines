cwlVersion: v1.0
class: CommandLineTool

requirements:
    DockerRequirement:
        dockerPull: "quay.io/biocontainers/picard:2.26.7--hdfd78af_0"

baseCommand: [picard, MarkDuplicates]

inputs:
  INPUT:
    type: File
    inputBinding:
      separate: false
      prefix: I=
      position: 1
  OUTPUT:
    type: string
    inputBinding:
      separate: false
      prefix: O=
      position: 2
  METRICS_FILE:    
    type: string?
    inputBinding:
      separate: false
      prefix: M=
      position: 3
  remove_duplicates:
    type: boolean?
    default: false
    inputBinding:
      prefix: --REMOVE_DUPLICATES
      position: 4
  remove_sequencing_duplicates:
    type: boolean?
    default: false
    inputBinding:
      prefix: --REMOVE_SEQUENCING_DUPLICATES
      position: 5

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.OUTPUT)
  output_report: 
    type: File
    outputBinding:
        glob: $(inputs.METRICS_FILE)
