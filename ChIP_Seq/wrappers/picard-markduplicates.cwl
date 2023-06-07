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
            position: 3
    OUTPUT:
        type: string?
        inputBinding:
            separate: false
            prefix: O=
            position: 4
    metrics:
        type: string?
        inputBinding:
            separate: false
            prefix: M=
            position: 5
    remove_duplicates:
        type: boolean
        default: false
        inputBinding:
            separate: false
            prefix: REMOVE_DUPLICATES=
            position: 6 
    validation_stringency:
        type: string?
        inputBinding:
            separate: false
            prefix: VALIDATION_STRINGENCY=
            position: 7

outputs:
    output:
        type: File
        outputBinding:
            glob: $(inputs.OUTPUT)

    output_metrics:
        type: File
        outputBinding:
            glob: $(inputs.metrics)
