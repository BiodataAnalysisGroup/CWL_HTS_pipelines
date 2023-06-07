cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: DockerRequirement
  dockerPull: "biowardrobe2/macs2:v2.1.1"
- class: InlineJavascriptRequirement

baseCommand: ["macs2","callpeak"]

arguments:  
  - valueFrom: $(inputs.treatment.nameroot)
    prefix: --name
    position: 1

inputs:
  call_summits:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --call-summits
    doc: 'If set, MACS will use a more sophisticated signal processing approach to
        find subpeak summits in each enriched peak region. DEFAULT: False '
  f:
    type: string?
    inputBinding:
      position: 1
      prefix: --format
    doc: '-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE}, --format
        {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE} Format of tag file,
        "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM"
        or "BOWTIE" or "BAMPE". The default AUTO option will let MACS decide which format
        the file is. Note that MACS can''t detect "BAMPE" or "BEDPE" format with "AUTO",
        and you have to implicitly specify the format for "BAMPE" and "BEDPE". DEFAULT:
        "AUTO".'
  cutoff_analysis:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --cutoff-analysis
    doc:  |
      While set, MACS2 will analyze number or total length of peaks that can be
      called by different p-value cutoff then output a summary table to help user
      decide a better cutoff. The table will be saved in NAME_cutoff_analysis.txt
      file. Note, minlen and maxgap may affect the results. WARNING: May take ~30
      folds longer time to finish.
  p:
    type: float?
    inputBinding:
      position: 1
      prefix: --pvalue
    doc:  |
      Pvalue cutoff for peak detection. DEFAULT: not set.  -q, and -p are mutually
      exclusive. If pvalue cutoff is  set, qvalue will not be calculated and reported
      as -1  in the final .xls file..
  p_file:
    type: File?
    inputBinding:
      position: 1
      prefix: --pvalue
      loadContents: True
      valueFrom: ${ return inputs.input.contents.split('\n')[0];}
    doc: |
      Pvalue cutoff for peak detection loaded from the first line of a file.
  q:
    type: float?
    inputBinding:
      position: 1
      prefix: --qvalue
    doc: |
      The qvalue (minimum FDR) cutoff to call significant regions. Default is 0.05
  q_file:
    type: File?
    inputBinding:
      position: 1
      prefix: --qvalue
      loadContents: True
      valueFrom: ${ return inputs.q_file.contents.split('\n')[0];}
    doc: |
      The qvalue (minimum FDR) cutoff to call significant regions loaded from the first line of a file. Default is 0.05
  bdg:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --bdg
    doc:  |
      Store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files
  treatment:
    type: File
    inputBinding:
      position: 2
      prefix: --treatment
    doc:  |
      Treatment sample file(s). If multiple files are given as -t A B C, then
      they will all be read and pooled together. IMPORTANT: the first sample will
      be used as the outputs basename.
  control:
    type: File?
    inputBinding:
      position: 2
      prefix: --control
    doc: |
      The control or mock data file. Please follow the same direction as for -t/--treatment.
  gsize:
    type: string?
    inputBinding:
      position: 2
      prefix: --gsize
    doc: |
      It's the mappable genome size or effective genome size which is defined as the genome size which can be sequenced.
  nomodel:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --nomodel
    doc: |
      While on, MACS will bypass building the shifting model.
  shift:
    type: int?
    inputBinding:
      position: 2
      prefix: --shift
    doc: |
      Note, this is NOT the legacy --shiftsize option which is replaced by --extsize! You can set an arbitrary shift in bp here.
  extsize:
    type: int?
    inputBinding:
      position: 2
      prefix: --extsize
    doc: |
      While '--nomodel' is set, MACS uses this parameter to extend reads in 5'->3' direction to fix-sized fragments
  broad:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --broad
    doc: |
      When this flag is on, MACS will try to composite broad regions in BED12
  broad-cutoff:
    type: float?
    inputBinding:
      position: 2
      prefix: --broad-cutoff
    doc: |
      Cutoff for broad region. This option is not available unless --broad is set. If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff. DEFAULT: 0.1

outputs:
  narrowPeak:
    type: File
    outputBinding:
      glob: "*_peaks.narrowPeak"
  xls:
    type: File
    outputBinding:
      glob: "*_peaks.xls"
  bed:
    type: File
    outputBinding:
      glob: "*_summits.bed"
  lambda:
    type: File?
    outputBinding:
      glob: "*_control_lambda.bdg"
  pileup:
    type: File?
    outputBinding:
      glob: "*_treat_pileup.bdg"
  broadPeak:
    type: File?
    outputBinding:
      glob: "*_peaks.broadPeak"
  gappedPeak:
    type: File?
    outputBinding:
      glob: "*_peaks.gappedPeak"
  model_r:
    type: File?
    outputBinding:
      glob: "*_model.r"
  cutoff:
    type: File?
    outputBinding:
      glob: "*_cutoff_analysis.txt"
