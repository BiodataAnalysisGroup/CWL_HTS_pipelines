cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "ubuntu:latest"

baseCommand: ["cat"]

inputs:
    headerFile: 
      type: File
      inputBinding: 
        position: 1
    countsFile:
      type: File
      inputBinding: 
        position: 2
    outputFile: 
      type: string

outputs:
    output:
        type: stdout

stdout: $(inputs.outputFile)
