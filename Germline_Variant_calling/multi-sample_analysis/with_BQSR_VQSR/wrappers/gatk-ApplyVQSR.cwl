cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, ApplyVQSR]

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
    variant:
        type: File
        secondaryFiles:
        - .tbi
        inputBinding:
            position: 2
            prefix: -V
    ts_filter_level: 
        type: float?
        inputBinding:
            position: 3
            prefix: --truth-sensitivity-filter-level
    tranches_file:
        type: File
        inputBinding:
            position: 4
            prefix: --tranches-file
    recal_file:
        type: File
        secondaryFiles: 
        - .tbi
        inputBinding: 
            position: 5
            prefix: --recal-file
    mode: 
        type: string
        inputBinding: 
            position: 6
            prefix: -mode
    intervals:
        type: File
        inputBinding:
            prefix: -L
            position: 7
            shellQuote: false
    exclude_intervals:
        type: File?
        inputBinding:
            prefix: -XL
            position: 8
            shellQuote: false
    output_name:
        type: string?
        default: "cohort.vqsr"

arguments:
- prefix: -O
  position: 9
  shellQuote: false
  valueFrom: $(inputs.output_name + "." + inputs.intervals.basename.split("-")[0] + "." + inputs.mode + ".vcf.gz")

outputs:
    output:
        type: File
        secondaryFiles: 
        - .tbi
        outputBinding:
            glob: $(inputs.output_name + "." + inputs.intervals.basename.split("-")[0] + "." + inputs.mode + ".vcf.gz")
