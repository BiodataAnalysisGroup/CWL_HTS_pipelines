cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/samtools:1.14--hb421002_0"

baseCommand: [samtools, view]
inputs:
  isbam:
    type: boolean
    default: false
    inputBinding:
      position: 2
      prefix: -b
    doc: |
      output in BAM format
  readswithoutbits:
    type: int?
    inputBinding:
      position: 1
      prefix: -F
    doc: |
      Do not output alignments with any bits set in FLAG present in the FLAG field
  readswithbits:
    type: int?
    inputBinding:
      position: 1
      prefix: -f
    doc: |
      only include reads with all bits set in INT set in FLAG [0]
  collapsecigar:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -B
    doc: |
      collapse the backward CIGAR operation
  readsingroup:
    type: string?
    inputBinding:
      position: 1
      prefix: -r
    doc: |
      only include reads in read group STR [null]
  uncompressed:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -u
    doc: |
      uncompressed BAM output (implies -b)
  readtagtostrip:
    type: string[]?
    inputBinding:
      position: 1
    doc: |
      read tag to strip (repeatable) [null]
  input:
    type: File
    inputBinding:
      position: 4
    doc: |
      Input bam file.
  readsquality:
    type: int?
    inputBinding:
      position: 1
      prefix: -q
    doc: |
      only include reads with mapping quality >= INT [0]
  cigar:
    type: int?
    inputBinding:
      position: 1
      prefix: -m
    doc: |
      only include reads with number of CIGAR operations
      consuming query sequence >= INT [0]
  iscram:
    type: boolean
    default: false
    inputBinding:
      position: 2
      prefix: -C
    doc: |
      output in CRAM format
  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: -@
    doc: |
      number of BAM compression threads [0]
  fastcompression:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: '-1'
    doc: |
      use fast BAM compression (implies -b)
  samheader:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -h
    doc: |
      include header in SAM output
  count:
    type: boolean
    default: false
    inputBinding:
      position: 1
      prefix: -c
    doc: |
      Instead of printing the alignments, only count them and print the total number
  randomseed:
    type: float?
    inputBinding:
      position: 1
      prefix: -s
    doc: |
      integer part sets seed of random number generator [0];
      rest sets fraction of templates to subsample [no subsampling]
  region:
    type: string?
    inputBinding:
      position: 5
    doc: |
      [region ...]
  readsinlibrary:
    type: string?
    inputBinding:
      position: 1
      prefix: -l
    doc: |
      only include reads in library STR [null]
  target_bed_file:
    type: File?
    inputBinding: 
      position: 1
      prefix: -L
    doc:  |
      Only output alignments overlapping the input BED FILE [null].
  output_name:
    type: string?
    inputBinding:
      position: 2
      prefix: -o
outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_name)