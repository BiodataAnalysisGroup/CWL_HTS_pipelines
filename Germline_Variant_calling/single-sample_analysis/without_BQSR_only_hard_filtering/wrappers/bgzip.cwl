cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: bioslimcontainers/tabix:1.7

baseCommand: [bgzip]

inputs: 
    compression:
        type: boolean
        default: true
        inputBinding:
            position: 1
            prefix: -c
    input: 
        type: File
        inputBinding:
            position: 2

stdout: $(inputs.input.basename + ".gz")

outputs:
    output:
        type: stdout
