cwlVersion: v1.0
class: CommandLineTool

requirements:
#   - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0
#   - class: InitialWorkDirRequirement
#     listing: 
#     - entry: "$({class: 'Directory', listing: []})"
#       entryname: $(inputs.genomicsdb)
#       writable: true

# $namespaces:
#     cwltool: "http://commonwl.org/cwltool#"
#     hints:
#     cwltool:LoadListingRequirement:
#         loadListing: shallow_listing

baseCommand: [gatk, GenotypeGVCFs]

# arguments:
# - shellQuote: false
#   prefix: -V
#   position: 1
#   valueFrom: |
#    ${
#     var genomicsdb = "gendb://" + inputs.genomicsdb.path;
#     return genomicsdb;
#    }

inputs:
    reference:
        type: File
        secondaryFiles:
        - .amb
        - .ann
        - .bwt
        - .pac
        - .sa
        - .fai
        - ^.dict
        inputBinding: 
            position: 1
            prefix: -R
    # genomicsdb: 
    #     type: Directory
    #     # inputBinding: 
    #     #     position: 2
    #     #     prefix: -V
    gvcf_input:
        type: File
        secondaryFiles:
        - .tbi
        inputBinding:
            position: 2
            prefix: -V
    vcf_output:
        type: string
        default: "cohort.genotyped.vcf.gz"
        inputBinding:
            position: 3
            prefix: -O

outputs:
    output:
        type: File
        outputBinding:
            glob: $(inputs.vcf_output)
        secondaryFiles: 
        - .tbi