cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "biocontainers/bcftools:v1.5_cv3"

baseCommand: [bcftools, norm]
inputs:
  input: 
    type: File
    inputBinding: 
        position: 100
  threads: 
    type: int?
    inputBinding: 
        prefix: --threads
        position: 1
  reference:
    type: File
    inputBinding: 
        position: 2
        prefix: -f
  multiallelics:
    type: string
    inputBinding: 
        position: 3
        prefix: --multiallelics
    doc: |
     split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). 
     An optional type string can follow which controls variant types which should be split or merged together: 
     If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; 
     if SNPs and indels should be merged into a single record, specify any.  
  output_type: 
    type: string?
    default: v
    inputBinding:
        separate: false
        position: 4
        prefix: -O

arguments: 
- prefix: -o
  valueFrom: $(inputs.input.basename.split(".vcf")[0] + ".norm.vcf")

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.input.basename.split(".vcf")[0] + ".norm.vcf")