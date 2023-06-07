cwlVersion: v1.0
class: CommandLineTool

requirements:
    InlineJavascriptRequirement: {}

hints:
- class: DockerRequirement
  dockerPull: ubuntu:latest

baseCommand: cp

inputs:
    input_qc_check:
        type: boolean?
        default: true
    input_trimming_check:
        type: boolean?
        default: true
    input_file:
        type: File?
        default: "null"
        inputBinding:
            position: 1
    run_type:
        type: string?
outputs: 
    renamed_file:
      type: File
      outputBinding:
        glob: $(inputs.run_type + inputs.input_file.basename)

arguments:
    - position: 2
      valueFrom: $(inputs.run_type + inputs.input_file.basename)