cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, VariantRecalibrator]

arguments:
- shellQuote: false
  position: 9
  prefix: -O
  valueFrom: $(inputs.output_name + "." + inputs.intervals.basename.split("-")[0] + ".vcf.gz")
- shellQuote: false
  position: 10
  prefix: --tranches-file
  valueFrom: $(inputs.output_name + "." + inputs.intervals.basename.split("-")[0] + ".tranches")
- shellQuote: false
  position: 1
  valueFrom: $("--resource" + inputs.arguments_1 + " " + inputs.resource_1.path)

inputs:
    truth_sensitivity_trance: 
        type: float[]?
    use_annotation:
        type: string[]
        inputBinding:
            prefix: -an
            itemSeparator: " -an "
            shellQuote: false
            position: 1
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
    mode:
        type: string
        inputBinding:
            position: 3
            prefix: -mode
            shellQuote: false
    arguments_1: 
        type: string
    resource_1:
        type: File
        secondaryFiles: 
        - .tbi
    arguments_2: 
        type: string?
    resource_2:
        type: File?
        secondaryFiles: 
        - .tbi
    arguments_3: 
        type: string?
    resource_3:
        type: File?
        secondaryFiles: 
        - .tbi
    arguments_4: 
        type: string?
    resource_4:
        type: File?
        secondaryFiles: 
        - .tbi
    arguments_5: 
        type: string?
    resource_5:
        type: File?
        secondaryFiles: 
        - .tbi
    intervals:
        type: File?
    exclude_intervals:
        type: File?
        inputBinding:
            prefix: -XL
            position: 6
            shellQuote: false
    trust_all_polymorphic:
        type: boolean?
        inputBinding:
            position: 7
            prefix: --trust-all-polymorphic
            shellQuote: false
    max_gaussians:
        type: int?
        inputBinding:
            position: 8
            prefix: --max-gaussians
    output_name:
        type: string?
        default: "cohort"

outputs:
    recal_table:
        type: File
        secondaryFiles: 
        - .tbi
        outputBinding: 
            glob: $(inputs.output_name + "." + inputs.intervals.basename.split("-")[0] + ".vcf.gz")
    tranches_file:
        type: File
        outputBinding: 
            glob: $(inputs.output_name + "." + inputs.intervals.basename.split("-")[0] + ".tranches")
