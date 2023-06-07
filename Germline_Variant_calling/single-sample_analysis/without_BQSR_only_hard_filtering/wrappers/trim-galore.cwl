cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7

inputs:
  command:
    type: string
    default: trim_galore
    inputBinding:
      position: 1
  cutadapt:
    type: string?
    default: /usr/local/bin/cutadapt
    inputBinding:
      position: 2
      prefix: --path_to_cutadapt
  paired:
    type: boolean?
    default: false
    inputBinding:
      prefix: --paired
      position: 3
  fq_files:
    type: File[]?
    inputBinding:
      position: 4
  results_path:
    type: string?
    default: .
    inputBinding:
      prefix: -o
      position: 5
  quality:
    type: int?
    default: 20
    inputBinding:
      prefix: -q
      position: 6
  strigency:
    type: int?
    default: 1
    inputBinding:
      prefix: --stringency
      position: 7
    doc:  Overlap with adapter sequence required to trim a sequence. 
          Defaults to a very stringent setting of 1, i.e. even a single 
          bp of overlapping sequence will be trimmed off from the 3' end 
          of any read.
  length:
    type: int?
    default: 20
    inputBinding:
      prefix: --length
      position: 8
  compression:
    type: boolean?
    default: false
    inputBinding:
      prefix: --gzip
      position: 9
  do_not_compress:
    type: boolean?
    default: true
    inputBinding:
      prefix: --dont_gzip
      position: 10
    doc:  Output files won't be compressed with GZIP. 
          This option overrides --gzip.
  adapter:
    type: string?
    inputBinding:
      prefix: -a
      position: 11
  adapter2:
    type: string?
    inputBinding:
      prefix: -a2
      position: 12
  error_rate:
    type: float?
    default: 0.1
    inputBinding:
      prefix: -e
      position: 13
  no_report:
    type: boolean?
    default: false
    inputBinding:
      prefix: --no_report_file
      position: 14
  no_warnings:
    type: boolean?
    default: false
    inputBinding:
      prefix: --suppress_warn
      position: 15
  clip_r1:
    type: int?
    inputBinding:
      prefix: --clip_R1
      position: 16
  clip_r2:
    type: int?
    inputBinding:
      prefix: --clip_R2
      position: 17
  three_prime_clip_R1:
    type: int?
    inputBinding:
      prefix: --three_prime_clip_R1
      position: 18
  three_prime_clip_R2:
    type: int?
    inputBinding:
      prefix: --three_prime_clip_R2
      position: 19
  max_n: 
    type: int?
    inputBinding:
      prefix: --max_n
      position: 20
    doc:  The total number of Ns (as integer) a read may contain before it will be removed altogether.
          In a paired-end setting, either read exceeding this limit will result in the entire
          pair being removed from the trimmed output files.
  trim_n: 
    type: int?
    inputBinding:
      prefix: --trim_n
      position: 21
    doc:  Removes Ns from either side of the read. 
          This option does currently not work in RRBS mode.
  trim_suffix:
    type: string?
    doc: |
     Trimmed file suffix to separate trimmed fastq from trimming report, in case starting files begin with ".fq". 
     Set to ".fq.gz" or ".fq" for compression equal to true or false, respectively.

outputs:
  trim_galore:
    type: File[]
    outputBinding:
      glob: $("*" + inputs.trim_suffix) #"*.fq*"
  trim_galore_report:
    type: File[]
    outputBinding:
      glob: "*.txt"