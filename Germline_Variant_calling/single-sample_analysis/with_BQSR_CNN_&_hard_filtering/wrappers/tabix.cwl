cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: bioslimcontainers/tabix:1.7
  - class: InitialWorkDirRequirement
    listing: [ $(inputs.input) ]

baseCommand: [tabix]

inputs: 
    input: 
        type: File
        inputBinding:
            valueFrom: $(self.basename)

outputs:
    output:
        type: File
        outputBinding: 
            glob: $(inputs.input.basename)
        secondaryFiles:
            - .tbi
