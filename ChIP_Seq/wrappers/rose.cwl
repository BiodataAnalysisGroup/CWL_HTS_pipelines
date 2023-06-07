cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: InitialWorkDirRequirement
  listing:
    - entry: "$({class: 'Directory', listing: []})"
      entryname: "gff"
      writable: true
    - entry: "$({class: 'Directory', listing: []})"
      entryname: "mappedGFF"
      writable: true
    - entry: $(inputs.annotation)
      writable: true
    - entry: $(inputs.ROSE_main)
    - entry: $(inputs.bamToGFF)
    - entry: $(inputs.bamToGFF_turbo)
    - entry: $(inputs.ROSE_callSuper)
    - entry: $(inputs.ROSE_geneMapper)
    - entry: $(inputs.ROSE_main_turbo)
- class: DockerRequirement
  dockerPull: biodataanalysisgroup/rose_main:v1.0

doc: |
  This CWL wrapper runs ROSE to get Super Enhancer regions.
  ROSE code from BitBucket should be cloned with git (git clone https://bitbucket.org/young_computation/rose.git) to the 'wrappers' directory before executing this CWL wrapper.
  Detailed documentation is available here: http://younglab.wi.mit.edu/super_enhancer_code.html

baseCommand: ['python', 'ROSE_main.py']

arguments:
  - prefix: -o
    valueFrom: "."

inputs:
  # ROSE scripts & annotation files
  annotation:
    type: Directory
    default: {class: Directory, location: ../wrappers/rose/annotation}
  ROSE_main:
    type: File
    default: {class: File, location: ../wrappers/rose/ROSE_main.py}
  bamToGFF:
    type: File
    default: {class: File, location: ../wrappers/rose/ROSE_bamToGFF.py}
  bamToGFF_turbo:
    type: File
    default: {class: File, location: ../wrappers/rose/ROSE_bamToGFF_turbo.py}
  ROSE_callSuper: 
    type: File
    default: {class: File, location: ../wrappers/rose/ROSE_callSuper.R}
  ROSE_geneMapper: 
    type: File
    default: {class: File, location: ../wrappers/rose/ROSE_geneMapper.py}
  ROSE_main_turbo: 
    type: File
    default: {class: File, location: ../wrappers/rose/ROSE_main_turbo.py}
  # Input files
  macs_peaks:
    type: File
    inputBinding:
      position: 1
      prefix: -i
    doc: |
     ".gff file (described above) of regions that were previously calculated to be enhancers. I.e. Med1-enriched regions identified using MACS."
  ranking_bam:
    type: File
    inputBinding:
      position: 2
      prefix: -r
    secondaryFiles:
    - .bai
    doc: |
     ".bam file to be used for ranking enhancers by density of this factor. I.e. Med1 ChIP-Seq reads."
  genome_build:
    type: string
    inputBinding:
      position: 3
      prefix: -g
  stitch_distance:
    type: int?
    inputBinding:
      position: 4
      prefix: -s
    doc: |
     "Μaximum distance between two regions that will be stitched together (Default: 12.5kb)"
  tss_distance:
    type: int?
    inputBinding:
      position: 9
      prefix: -t
    doc: |
     "Εxclude regions contained within +/- this distance from TSS in order to account for promoter biases (Default: 0; recommended if used: 2500). If this value is 0, will not look for a gene file."
  control_bam:
    type: File?
    secondaryFiles: 
    - .bai
    inputBinding:
      position: 10
      prefix: -c
    doc: |
     ".bam file to be used as a control. Subtracted from the density of the RANKING_BAM. I.e. Whole cell extract reads."

outputs:
  gff_dir_outputs:
    type: File[]?
    outputBinding: 
      glob: "gff/*.gff"
  mappedGFF_dir_outputs:
    type: File[]?
    outputBinding: 
      glob: "mappedGFF/*.gff"
  STITCHED_ENHANCER_REGION_MAP:
    type: File
    outputBinding: 
      glob: "*_ENHANCER_REGION_MAP.txt"
  AllEnhancers_table: 
    type: File
    outputBinding: 
      glob: "*_AllEnhancers.table.txt"
  SuperEnhancers_table: 
    type: File
    outputBinding: 
      glob: "*_SuperEnhancers.table.txt"
  Enhancers_withSuper: 
    type: File
    outputBinding: 
      glob: "*_Enhancers_withSuper.bed"
  Plot_points:
    type: File
    outputBinding: 
      glob: "*_Plot_points.png"