cwlVersion: v1.0
class: CommandLineTool

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastqc:v0.11.5

inputs:
    command: 
        type: string
        default: /opt/FastQC/fastqc
        inputBinding:
            position: 1
    input_files:
        type: File[]
        inputBinding:
            position: 2
    fastqc_output_dir:
        type: string?
        default: .
        inputBinding:
            prefix: -o
            separate: true
            position: 3
outputs: 
    zipped_file:
      type:
        - File[]
      outputBinding:
        glob: '*.zip'
    html_file:
      type:
        - File[]
      outputBinding:
        glob: '*.html'
