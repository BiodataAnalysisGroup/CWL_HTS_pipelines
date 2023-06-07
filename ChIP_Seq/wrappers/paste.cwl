cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:latest"
- class: InlineJavascriptRequirement

baseCommand: [paste]

inputs:
    input_files:
        type: File[]?
        inputBinding:
            position: 1
    coordinates: 
        type: File?
        inputBinding: 
            position: 2
    counts: 
        type: File?
        inputBinding: 
            position: 3
    output_name:
        type: string

outputs:
    output:
        type: stdout

stdout: $(inputs.output_name)
