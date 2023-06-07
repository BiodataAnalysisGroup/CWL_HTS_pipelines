cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.14--hb421002_0"

baseCommand: [samtools, markdup]

inputs:
    threads:
        type: int?
        inputBinding:
            position: 1
            prefix: -@
    remove_duplicates: 
        type: boolean?
        inputBinding: 
            position: 2
            prefix: -r
        doc: |
            Remove duplicate reads
    input:
        type: File
        inputBinding:
            position: 100
        doc: |
            Input bam file
    output_name:
        type: string
        inputBinding:
            position: 101
outputs:
    output:
        type: File
        outputBinding:
            glob: "*.bam"