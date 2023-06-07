cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: broadinstitute/gatk:4.3.0.0
- class: InitialWorkDirRequirement
  listing: 
  - entry: "$({class: 'Directory', listing: []})"
    entryname: $(inputs.outputdir)
    writable: true

baseCommand: [gatk, SplitIntervals]

inputs:
    reference:
        type: File
        inputBinding:
            position: 1
            prefix: -R
        secondaryFiles:
        - .amb
        - .ann
        - .bwt
        - .pac
        - .sa
        - .fai
        - ^.dict
    include_intervalList:
        type: File?
        inputBinding: 
            position: 2
            prefix: -L
    exclude_intervalList:
        type: File?
        inputBinding: 
            position: 3
            prefix: -XL
    scatter_count:
        type: int?
        default: 1
        inputBinding: 
            position: 4
            prefix: --scatter-count
        label: Desired split for variant calling
    subdivision_mode:
      type: string?
      default: BALANCING_WITHOUT_INTERVAL_SUBDIVISION
      inputBinding:
        position: 5
        prefix: --subdivision-mode
    outputdir: 
        type: string?
        default: "interval_files_folder"

arguments: 
- prefix: "-O"
  valueFrom: $(inputs.outputdir)

outputs:
  # intervalfilesdir:
  #   type: Directory
  #   outputBinding: 
  #     glob: $(inputs.outputdir)
  intervalfiles:
    type: File[]
    label: Scatter intervals files
    outputBinding:
      glob: "$(inputs.outputdir)/*.interval_list"
