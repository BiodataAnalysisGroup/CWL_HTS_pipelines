cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:latest"
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

baseCommand: [awk]

arguments:
    - '$1~/^chr/ {print $0}'

inputs:
    input: 
        type: File
        inputBinding:
            position: 2

outputs:
    output:
        type: stdout

stdout: $(inputs.input.nameroot + "_processed.bed")
