cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "staphb/bedtools:2.30.0"
- class: InlineJavascriptRequirement

baseCommand: [bedtools, merge]
inputs:
    input:
        type: File
        inputBinding:
            prefix: -i
            position: 1
    output_name:
        type: string

outputs:
    output:
        stdout

stdout: $(inputs.output_name)
