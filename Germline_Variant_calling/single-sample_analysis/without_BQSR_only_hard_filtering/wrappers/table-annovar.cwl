cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: bioinfochrustrasbourg/annovar:2018Apr16

baseCommand: table_annovar.pl

inputs:
  query_file:
    type: File
    inputBinding:
      position: 1
  database_location:
    type: Directory
    inputBinding:
      position: 2
  build_over:
    type: string
    default: hg19
    inputBinding:
      prefix: --buildver
      position: 3
  remove:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --remove
  protocol:
    type: string
    inputBinding:
      prefix: --protocol
      position: 6
  operation:
    type: string
    inputBinding:
      prefix: --operation
      position: 7
  na_string:
    type: string?
    default: .
    inputBinding:
      prefix: --nastring
      position: 8
  vcfinput:
    type: boolean
    default: true
    inputBinding:
      position: 9
      prefix: --vcfinput
  otherinfo:
    type: boolean?
    inputBinding:
      position: 10
      prefix: --otherinfo
  convert_arg:
    type: string?
    inputBinding:
      position: 11
      prefix: --convertarg
  csvout:
    type: float?
    inputBinding:
      position: 12
      prefix: -csvout
  threads:
    type: int?
    inputBinding:
      position: 13
      prefix: --thread
  verbose:
    type: float?
    inputBinding:
      position: 14
      prefix: --verbose
  xreffile:
    type: File?
    inputBinding:
      position: 15
      prefix: --xreffile
  arguments_str:
    type: string?
    inputBinding: 
      position: 16
      prefix: --argument
  dbtype: 
    type: string?
    inputBinding: 
      position: 17
      prefix: --dbtype
  gff3dbfile:
    type: string?
    inputBinding:
      position: 18
      prefix: --gff3dbfile
  bedfile:
    type: string?
    inputBinding:
      position: 19
      prefix: --bedfile
  vcfdbfile:
    type: string?
    inputBinding:
      position: 20
      prefix: --vcfdbfile
  # output_name:
  #   type: string
  #   inputBinding:
  #     position: 4
  #     prefix: --out

arguments:
- prefix: --out
  valueFrom: $("annotated_" + inputs.query_file.basename.split(".vcf")[0])

outputs:
  multianno_vcf:
    type: File
    outputBinding:
      glob: "*.vcf"
  multianno_txt:
    type: File
    outputBinding:
      glob: "*.txt"
  avinput:
    type: File
    outputBinding:
      glob: "*.avinput"