cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - '$({class: "Directory", basename: inputs.outputdir, listing: []})'
  DockerRequirement:
    dockerPull: biodataanalysisgroup/chipqc-diffbind-rscripts:v1.0

baseCommand: ["Rscript", "--vanilla"]

inputs:
    rscript: 
        type: string
        default: /home/analysis/ChIPQC.R
        inputBinding:
            position: 1
    treatmentBAM: 
        type: File[]
        secondaryFiles: 
            - .bai
        inputBinding: 
            position: 2
            prefix: --treatment-bam
    controlBAM: 
        type: File[]?
        secondaryFiles: 
            - .bai
        inputBinding: 
            position: 3
            prefix: --control-bam
    peaks: 
        type: File[]
        inputBinding:
            position: 4
            prefix: --peak-files
    peakCaller: 
        type: string
        default: bed
        inputBinding: 
            position: 5
            prefix: --peak-caller
    metadata:
        type: File
        inputBinding:
            position: 6
            prefix: --metadata
    consensus:
        type: boolean?
        inputBinding: 
            position: 7
            prefix: --consensus
    bCount:
        type: boolean?
        inputBinding: 
            position: 8
            prefix: --bCount
    summits:
        type: int?
        default: 200
        inputBinding:
            position: 9
            prefix: --summits
    annotation:
        type: string
        default: "hg38"
        inputBinding:
            position: 10
            prefix: --annotation
    chromosomes_str:
        type: string?
        inputBinding:
            position: 11
            prefix: --chromosomes
        doc: |
         Specify one chromosome name for sampling (e.g., "chr18") or set value to "all" to use all chromosomes. Incompatible with "chromosomes_int" argument.
    chromosomes_int:
        type: int[]?
        inputBinding:
            position: 12
            prefix: --chromosomes
        doc: |
         Specify integer index for chromosomes to be used for sampling (e.g., [1, 2, 10]). Incompatible with "chromosomes_str" argument.
    blacklist:
        type: File?
        inputBinding:
            position: 12
            prefix: --blacklist
    mapQCth:
        type: int?
        default: 15
        inputBinding:
            position: 13
            prefix: --mapQCth
    profileWin:
        type: int?
        default: 400
        inputBinding:
            position: 14
            prefix: --profileWin
    fragmentLength:
        type: int?
        default: 125
        inputBinding:
            position: 15
            prefix: --fragmentLength
    bParallel:
        type: boolean?
        default: true
        inputBinding:
            position: 16
            prefix: --bParallel
    facetBy:
        type: string[]?
        inputBinding:
            position: 17
            prefix: --facetBy
    outputdir:
        type: string
        default: "ChIPQC_HTML_report"
        inputBinding: 
            position: 18
            prefix: --reportFolder
    chipqc_experiment:
        type: string
        inputBinding: 
            position: 19
            prefix: --chipqc-experiment
    chipqc_report:
        type: string
        inputBinding: 
            position: 20
            prefix: --chipqc-report

outputs:
    ChIPQCexperiment:
        type: File
        outputBinding:
            glob: $(inputs.chipqc_experiment)
    outdir:
        type: Directory?
        outputBinding:
            glob: $(inputs.outputdir)
    ChIPQCreport:
        type: File?
        outputBinding:
            glob: $(inputs.chipqc_report)
