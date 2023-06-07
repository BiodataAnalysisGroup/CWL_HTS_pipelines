cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: 
    - entry: "$({class: 'Directory', listing: []})"
      entryname: $(inputs.outputdir)
      writable: true

hints:
  DockerRequirement:
    dockerPull: ubuntu:latest

arguments:
    - position: 1
      valueFrom: $(inputs.command)
    - position: 2
      valueFrom: $("-t")
    - position: 3
      valueFrom: $(inputs.outputdir)

inputs:
    input_qc_check:
        type: boolean?
        default: true
    input_trimming_check:
        type: boolean?
        default: true
    command:
      type: string?
      default: cp
    input_files:
        # type: File?
        # default: "null"
        # type: 
        # - "null"
        # - type: array
        #   items: File?
        type: 
        - type: array
          items: ["null", File]
        inputBinding:
            position: 10
    outputdir:
        type: string?

outputs: 
    output_dir:
        type: Directory?
        outputBinding:
            glob: $(inputs.outputdir)
