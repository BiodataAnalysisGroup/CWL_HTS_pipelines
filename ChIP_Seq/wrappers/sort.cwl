cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:latest"
- class: InlineJavascriptRequirement

baseCommand: [sort]

inputs:
    unique:
        type: boolean?
        inputBinding:
            position: 1
            prefix: -u
    key1: 
        type: string
        inputBinding: 
            position: 2
            prefix: -k
            separate: false
            shellQuote: false
    key2: 
        type: string?
        inputBinding: 
            position: 3
            prefix: -k
            separate: false
            shellQuote: false
    input: 
        type: File
        inputBinding:
            position: 4
    output_name:
        type: string

outputs:
    output:
        type: stdout

stdout: $(inputs.output_name)
