cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0"
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement

arguments:
- cat
- $(inputs.input)
- shellQuote: false
  valueFrom: '|'
- awk
- '{print $1 "," $2 ","  $3}'
- shellQuote: false
  valueFrom: '|'
- csvtk
- csv2tab

inputs:
    input: 
      type: File[]
    consensus_bed:
      type: string
outputs:
    output:
        type: stdout

stdout: $(inputs.consensus_bed)
