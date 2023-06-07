cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:latest"
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

doc: |
 For compatibility with ROSE, the resulting .gff must have the following columns:
    1: chromosome (chr#)
    2: unique ID for each constituent enhancer region
    4: start of constituent
    5: end of constituent
    7: strand (+,-,.)
    9: unique ID for each constituent enhancer region

baseCommand: [awk]

arguments:
    - position: 1
      valueFrom: '{print $1,$4,"MACS",$2,$3,$10,$6,".",$4}'
    - position: 2
      shellQuote: false
      valueFrom: "OFS='\t'"      

inputs:
    input: 
        type: File
        inputBinding:
            position: 3

outputs:
    output:
        type: stdout

stdout: $(inputs.input.nameroot + "_rose.gff")
