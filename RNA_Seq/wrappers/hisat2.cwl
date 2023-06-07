cwlVersion: v1.0
class: CommandLineTool

baseCommand: hisat2

hints:
- class: DockerRequirement
  dockerPull: greatfireball/hisat2

requirements:
- class: InlineJavascriptRequirement

arguments:
- "-x"
- $(inputs.idx_directory.path)/$(inputs.idx_basename)

inputs:
  num_of_threads:
    type: int?
    default: 16
    inputBinding:
      prefix: -p
      position: 1
    doc:  The -p option causes HISAT2 to launch a specified number of parallel search threads. 
          Each thread runs on a different processor/core and all threads find alignments in parallel, 
          increasing alignment throughput by approximately a multiple of the number of threads (though in practice, speedup is somewhat worse than linear).
  alignments_tailored_trans_assemb:
    type: boolean
    default: true
    inputBinding:
      prefix: --dta
      position: 2
    doc:  Report alignments tailored for transcript assemblers including StringTie. 
          With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. 
          This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage.
  idx_directory:
    type: Directory
  idx_basename:
    label: "Basename of the hisat2 index files"
    doc: "Basename of the hisat2 index files, not including extensions like .1.ht2"
    type: string
  files_with_first_mates:
    type: File?
    inputBinding:
      prefix: "-1"
      position: 3
    doc:  Comma-separated list of files containing mate 1s (filename usually includes _1), e.g. -1 flyA_1.fq,flyB_1.fq. 
          Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m2>. 
          Reads may be a mix of different lengths. If - is specified, hisat2 will read the mate 1s from the “standard in” or “stdin” filehandle.
  files_with_second_mates:
    type: File?
    inputBinding:
      prefix: "-2"
      position: 4
    doc:  Comma-separated list of files containing mate 2s (filename usually includes _2), e.g. -2 flyA_2.fq,flyB_2.fq. 
          Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m1>. 
          Reads may be a mix of different lengths. If - is specified, hisat2 will read the mate 2s from the “standard in” or “stdin” filehandle.
  files_with_unpaired_reads:
    type: File?
    inputBinding:
      prefix: -U
      position: 5
    doc:  Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. 
          Reads may be a mix of different lengths. If - is specified, hisat2 gets the reads from the “standard in” or “stdin” filehandle.
  SAM_output:
    type: string
    inputBinding:
      prefix: -S
      position: 6
  fastq_format:
    type: boolean?
    default: true
    inputBinding: 
      prefix: -q
      position: 7
  known_splicesite_infile:
    type: File?
    inputBinding:
      prefix: --known-splicesite-infile
      position: 8
  unpaired_reads_unaligned:
    type: string?
    inputBinding:
      prefix: --un-gz
      position: 9
  paired_reads_not_concordant:
    type: string?
    inputBinding:
      prefix: --un-conc-gz
      position: 10
  remove_chrname:
    type: boolean?
    default: false
    inputBinding:
      prefix: --remove-chrname
      position: 11
  add_chrname:
    type: boolean?
    default: false
    inputBinding:
      prefix: --add-chrname
      position: 12

stderr: $(inputs.SAM_output.split(".sam")[0].concat("_hisat2.report.txt"))

outputs:
  output:
    type: File
    outputBinding:
      glob: "*.sam"
  output_stderr:
    type: stderr
