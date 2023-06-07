cwlVersion: v1.0
class: CommandLineTool

requirements:
    - class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: ubuntu:latest

baseCommand: cp

inputs:
    input_file:
        type: File
        inputBinding:
            position: 1
    run_type:
        type: string
outputs: 
    renamed_file:
      type: File
      outputBinding:
        glob: $(inputs.run_type + inputs.input_file.basename)

arguments:
    - position: 2
      valueFrom: $(inputs.run_type + inputs.input_file.basename)