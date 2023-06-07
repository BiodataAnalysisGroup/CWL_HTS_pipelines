cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:latest"
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

baseCommand: [awk]

inputs:
    print: 
        type: string
        inputBinding:
            position: 1
    ofs:
        type: string
        default: '\t'
        inputBinding:
            position: 2
            prefix: OFS=
            separate: false
    input_file: 
        type: File
        inputBinding:
            position: 3
    second_file: 
        type: File?
        inputBinding:
            position: 4
    split_string:
        type: string
    output_suffix:
        type: string

outputs:
    output:
        type: stdout

stdout: $(inputs.input_file.basename.split(inputs.split_string)[0].concat(inputs.output_suffix))
