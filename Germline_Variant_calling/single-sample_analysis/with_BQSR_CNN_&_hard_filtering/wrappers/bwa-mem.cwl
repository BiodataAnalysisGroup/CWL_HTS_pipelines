cwlVersion: v1.0
class: CommandLineTool
doc: "[bwa](https://github.com/lh3/bwa)"

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"

baseCommand: bwa
arguments: ["mem"]

inputs:
  file_split:
    type: string?
    default: "_R"
  sec_shorter_split_hits:
    type: boolean
    default: true
    inputBinding:
      prefix: -M
      position: 1
    doc: "Mark shorter split hits as secondary (for Picard compatibility)."
  num_threads:
    type: int
    default: 16
    inputBinding:
      prefix: -t
      position: 2
    doc: "Number of threads"
  min_seed_length:
    type: int?
    inputBinding:
      prefix: -k
      position: 3
    doc: |
     "Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [19]"
  band_width:
    type: int?
    inputBinding:
      prefix: -w
      position: 4
    doc: |
     "Band width. Essentially, gaps longer than INT will not be found. 
     Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. [100]"
  dropoff:
    type: int?
    inputBinding:
      prefix: -d
      position: 5
    doc: |
     "Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, 
     where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff 
     except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, 
     but also reduces poor alignments inside a long good alignment. [100]"
  reseeding: 
    type: float?
    inputBinding:
      prefix: -r
      position: 6
    doc: |
     "Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. 
     Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [1.5]"
  discard_mem:
    type: int?
    inputBinding:
      prefix: -c
      position: 7
    doc: |
     "Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [10000]"
  P: 
    type: float?
    inputBinding:
      prefix: -P
      position: 8
    doc: |
     "In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair."
  matching_score: 
    type: int?
    inputBinding: 
      prefix: -A
      position: 9
    doc: "Matching score. [1]"
  mismatch_penalty:
    type: int?
    inputBinding:
      prefix: -B
      position: 10
    doc: "Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4]"
  gap_open_penalty:
    type: int?
    inputBinding:
      prefix: -O
      position: 11
    doc: "Gap open penalty. [6]"
  gap_extension_penalty:
    type: int?
    inputBinding:
      prefix: -E
      position: 12
    doc: "Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [1]"
  clipping_penalty: 
    type: int?
    inputBinding:
      prefix: -L
      position: 13
    doc: |
     "Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. 
     If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, 
     the SAM AS tag reports the best SW score; clipping penalty is not deducted. [5]"
  unpaired_read_pair_penalty:
    type: int?
    inputBinding:
      prefix: -U
      position: 14
    doc: |
     "Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. 
     It compares these two scores to determine whether we should force pairing. [9]"
  p: 
    type: float?
    inputBinding:
      prefix: -p
      position: 15
    doc: |
     "Assume the first input query file is interleaved paired-end FASTA/Q. See the command description for details."
  R: 
    type: string?
    inputBinding: 
      prefix: -R
      position: 16
    doc: |
     "Complete read group header line. ’\t’ can be used in STR and will be converted to a TAB in the output SAM. The read group ID will be attached to every read in the output. An example is ’@RG\tID:foo\tSM:bar’. [null]"
  T: 
    type: int?
    inputBinding:
      prefix: -T
      position: 17
    doc: "Don’t output alignment with score lower than INT. This option only affects output. [30]"
  a: 
    type: float? 
    inputBinding:
      prefix: -a
      position: 18
    doc: "Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments."
  C: 
    type: float? 
    inputBinding:
      prefix: -C
      position: 19
    doc: |
     "Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output. 
     Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). 
     Malformated comments lead to incorrect SAM output."
  hard_clipping: 
    type: float?
    inputBinding:
      prefix: -H
      position: 20
    doc: |
     "Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences."
  verbose:
    type: int?
    inputBinding: 
      prefix: -v
      position: 21
    doc: |
     "Control the verbose level of the output. This option has not been fully supported throughout BWA. Ideally, a value 0 for disabling all the output to stderr; 
     1 for outputting errors only; 2 for warnings and errors; 3 for all normal messages; 4 or higher for debugging. 
     When this option takes value 4, the output is not SAM. [3]"
  ref_genome:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
    inputBinding:
      position: 100
  trimmed_fq_read1:
    type: File
    inputBinding:
      position: 101
  trimmed_fq_read2:
    type: File?
    inputBinding:
      position: 102

stdout: $(inputs.trimmed_fq_read1.basename.split(inputs.file_split)[0].concat(".sam"))

outputs:
  output:
    type: stdout