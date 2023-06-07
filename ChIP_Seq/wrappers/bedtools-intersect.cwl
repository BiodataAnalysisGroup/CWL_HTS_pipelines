cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "staphb/bedtools:2.30.0"
- class: InlineJavascriptRequirement

baseCommand: [bedtools, intersect]

inputs:
    no_overlap:
        type: boolean?
        default: false
        inputBinding:
            prefix: -v
            position: 1
    wa:
        type: boolean?
        default: false
        inputBinding:
            prefix: -wa
            position: 2
    feature_a:
        type: File
        inputBinding:
            prefix: -a
            position: 3
    feature_b:
        type: 
        - File
        - File[]
        inputBinding:
            prefix: -b
            position: 4   
    outputFile:
        type: string
outputs:
    output:
        stdout

stdout: $(inputs.outputFile)
