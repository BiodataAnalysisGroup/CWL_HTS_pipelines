cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:latest"
- class: InlineJavascriptRequirement

baseCommand: [awk]

arguments:
- 'OFS="\t" {print $1,$2,$3}'

inputs:
    input_file:
        type: File
        inputBinding:
            position: 2

outputs:
    output:
        type: stdout

stdout: $(inputs.input_file.basename.split("_markdup_counts.bed")[0].concat(".coord"))
