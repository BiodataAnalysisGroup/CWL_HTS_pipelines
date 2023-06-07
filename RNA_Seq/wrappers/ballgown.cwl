cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

hints:
- class: DockerRequirement
  dockerPull: biodataanalysisgroup/ballgown-rscript:v1.0

baseCommand: ["Rscript", "--vanilla"]

inputs:
    stringtie_dirs:
        type: Directory[]
        inputBinding:
            position: 100
            prefix: --stringtie-output
    rscript: 
        type: string
        default: /home/analysis/ballgown.R
        inputBinding:
            position: 1
    phenotype_file:
        type: File
        inputBinding:
            position: 2
            prefix: --metadata
        doc: TSV/CSV file containing the sample names, the phenotype(s) of interest and the potential confounders
    phenotype:
        type: string
        inputBinding:
            position: 3
            prefix: --phenotype
        doc: Name of the column with phenotype of interest
    samples:
        type: string
        inputBinding:
            position: 4
            prefix: --samples
        doc: Name of the column with sample names/IDs
    timecourse:
        type: boolean?
        default: false
        inputBinding:
            position: 5
            prefix: --timecourse
    feature: 
        type: string
        default: transcript
        inputBinding:
            position: 6
            prefix: --feature
    measure: 
        type: string
        default: FPKM
        inputBinding: 
            position: 7
            prefix: --measure
    confounders:
        type: string?
        inputBinding:
            position: 8
            prefix: --confounders
        doc: |
         "optional comma-separated string (split into a vector of strings within Rscript) representing the names of potential confounders. 
         Must correspond to names of columns of pData(gown), which here are extracted from the phenotype_file."
    custom_model:
        type: boolean?
        default: false
        inputBinding:
            position: 9
            prefix: --custom-model
        doc: "(Optional) provide the design matrices for the models to be compared"
    mod: 
        type: string?
        inputBinding: 
            position: 10
            prefix: --mod
        doc: |
         "optional comma-separated string (split into a vector of strings within Rscript) representing the names of variables 
         (including variable(s) to be tested in custom models)"
    mod0: 
        type: string?
        inputBinding: 
            position: 11
            prefix: --mod0
        doc: |
         "optional comma-separated string (split into a vector of strings within Rscript) representing the names of variables 
         (excluding variable(s) to be tested in custom models)"

outputs:
    ballgown_de_results:
        type: File
        outputBinding:
            glob: "ballgown_DE_results.tsv"
    ballgown_object:
        type: File
        outputBinding:
            glob: "ballgown.rda"
    ballgown_de_custom_model:
        type: File?
        outputBinding:
            glob: "ballgown_CM_DE_results.tsv"
