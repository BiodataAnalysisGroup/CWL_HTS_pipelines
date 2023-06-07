cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: genomicpariscentre/featurecounts

baseCommand: featureCounts

inputs:
  annotation_file:
    type: File
    inputBinding:
      prefix: -a
      position: 1
    doc:  |
      Name of an annotation file. GTF/GFF format by default. 
      See  -F  option  for  more format  information.  
      Inbuilt annotations (SAF format) is available in 'annotation' directory of the package.
  number_of_threads:
    type: int?
    default: 16
    inputBinding:
      prefix: -T
      position: 2
  format_type:
    type: string?
    default: "GTF"
    inputBinding: 
      prefix: -F
      position: 3
  feature_type:
    type: string?
    default: "exon"
    inputBinding:
      prefix: -t
      position: 4
    doc:  |
      Specify feature type in GTF annotation. 'exon' by default. 
      Features used for read counting will be extracted from annotation using the provided value.
  attribute_type: 
    type: string?
    default: "gene_id"
    inputBinding:
      prefix: -g
      position: 5
    doc:  |
      Specify  attribute type in GTF annotation. 'gene_id' by default. 
      Meta-features used for read counting will be extracted from annotation using the provided value.
  chromosome_name_alias:
    type: File?
    inputBinding:
      prefix: -A
      position: 6
    doc:  |
      Provide a chromosome name alias file to match chr names in annotation with those in
      the  reads.  This should be a twocolumn comma-delimited text file. Its first column
      should include chr names in the annotation and its second column should include chr
      names  in  the  reads.  Chr  names  are  case sensitive. No column header should be
      included in the file.
  feature_level_counting:
    type: boolean?
    inputBinding:
      prefix: -f
      position: 7
    doc:  |
      Perform read counting at feature level (eg. counting reads for  exons  rather  than
      genes).
  read_meta_feature_overlap:
    type: boolean?
    inputBinding:
      prefix: -O
      position: 8
    doc:  |
      Assign  reads  to  all  their  overlapping  meta-features  (or  features  if  -f is
      specified).
  minOverlap: 
    type: int?
    default: 1
    inputBinding:
      prefix: --minOverlap
      position: 9
    doc:  |
      Minimum number of overlapping bases in a read that is required for read assignment.
      1 by default. Number of overlapping bases is counted from both reads if paired end.
      If a negative value is provided, then a gap of up to specified size will be allowed
      between read and the feature that the read is assigned to.
  fracOverlap: 
    type: float?
    default: 0
    inputBinding:
      prefix: --fracOverlap
      position: 10
    doc:  |
      Minimum fraction of overlapping bases in a read that is
      required  for  read  assignment.  Value should be within range [0,1]. 0 by default.
      Number of overlapping bases is counted from both reads if  paired  end.  Both  this
      option and '--minOverlap' option need to be satisfied for read assignment.
  largestOverlap: 
    type: boolean?
    inputBinding:
      prefix: --largestOverlap
      position: 12
    doc:  |
      Assign reads to a meta-feature/feature that has the largest number  of  overlapping
      bases.
  readExtension5: 
    type: boolean?
    inputBinding:
      prefix: --readExtension5
      position: 13
    doc:  |
      Reads are extended upstream by <int> bases from their
      5' end.
  readExtension3: 
    type: boolean?
    inputBinding:
      prefix: --readExtension3
      position: 14
    doc:  |
      Reads are extended upstream by <int> bases from their
      3' end.
  read2pos:
    type: boolean?
    inputBinding:
      prefix: --read2pos
      position: 15
    doc:  |
      Reduce reads to their 5' most base or 3' most base. Read counting is then performed
      based on the single base the read is reduced to.
  multi_mapping_reads:
    type: boolean?
    inputBinding:
      prefix: -M
      position: 16
    doc:  |
      Multi-mapping reads will also be counted. For a multimapping read, all its reported
      alignments  will  be  counted.  The  'NH'  tag  in  BAM/SAM input is used to detect
      multi-mapping reads.
  fraction:
    type: boolean?
    inputBinding:
      position: 17
      prefix: --fraction
    doc:  |
      Assign fractional counts to features. This option must be used together  with  '-M'
      or  '-O'  or  both.  When  '-M'  is  specified,  each  reported  alignment  from  a
      multi-mapping read (identified via 'NH' tag) will carry a fractional count of  1/x,
      instead of 1 (one), where x is the total number of alignments reported for the same
      read. When '-O' is specified, each overlapping feature will  receive  a  fractional
      count  of  1/y,  where y is the total number of features overlapping with the read.
      When both '-M' and '-O' are specified, each alignment will carry a fractional count
      of 1/(x*y).
  quality:
    type: int?
    default: 0
    inputBinding:
      position: 18
      prefix: -Q
    doc:  |
      The  minimum  mapping quality score a read must satisfy in order to be counted. For
      paired-end reads, at least one end should satisfy this criteria. 0 by default.
  splitOnly:
    type: boolean?
    inputBinding:
      position: 19
      prefix: --splitOnly
    doc:  
      Count split alignments only (ie. alignments with CIGAR string containing  'N').  An
      example of split alignments is exon-spanning reads in RNA-seq data.
  nonSplitOnly:
    type: boolean?
    inputBinding:
      prefix: --nonSplitOnly
      position: 20
    doc:  |
      If  specified,  only non-split alignments (CIGAR strings do not contain letter 'N')
      will be counted. All the other alignments will be ignored.
  primary:
    type: boolean?
    inputBinding:
      prefix: --primary
      position: 21
    doc:  |
      Count primary alignments only. Primary alignments are identified using bit 0x100 in
      SAM/BAM FLAG field.
  ignoreDup:
    type: boolean?
    inputBinding: 
      prefix: --ignoreDup
      position: 22
    doc:  |
      Ignore  duplicate  reads in read counting. Duplicate reads are identified using bit
      Ox400 in BAM/SAM FLAG field. The whole read pair is ignored if one of the reads  is
      a duplicate read for paired end data.
  strandness:
    type: int?
    default: 0
    inputBinding:
      prefix: -s
      position: 23
    doc:  |
      Perform  strand-specific  read  counting.  Acceptable  values:  0  (unstranded),  1
      (stranded) and 2 (reversely stranded).  0 by default.
  exon_junction_support:
    type: boolean?
    inputBinding: 
      prefix: -J
      position: 24
    doc:  |
      Count  number  of  reads  supporting  each  exon-exon  junction.   Junctions   were
      identified  from  those  exon-spanning  reads in the input (containing 'N' in CIGAR
      string). Counting results are saved to a file named '<output_file>.jcounts'
  genome_fasta:
    type: File?
    inputBinding:
      prefix: -G
      position: 25
    doc:  |
      Provide the name of a FASTA-format file that contains the reference sequences  used
      in  read  mapping  that produced the provided SAM/BAM files. This optional argument
      can be used with '-J' option to improve read counting for junctions.
  count_fragments:
    type: boolean?
    inputBinding:
      prefix: -p
      position: 26
    doc:  |
      If specified, fragments (or templates) will  be  counted  instead  of  reads.  This
      option is only applicable for paired-end reads.
  B: 
    type: boolean?
    inputBinding: 
      prefix: -B
      position: 27
    doc:  |
      Only count read pairs that have both ends aligned.
  P: 
    type: boolean?
    inputBinding:
      prefix: -P
      position: 28
    doc:  |
      Check  validity  of  paired-end distance when counting read pairs. Use -d and -D to
      set thresholds.
  min_frag_length:
    type: int?
    default: 50
    inputBinding:
      position: 29
      prefix: -d
    doc:  |
      Minimum fragment/template length, 50 by default.
  max_frag_length:
    type: int?
    default: 600
    inputBinding:
      position: 30
      prefix: -D
    doc:  |
      Maximum fragment/template length, 50 by default.
  C: 
    type: boolean?
    inputBinding:
      position: 31
      prefix: -C
    doc:  |
      Do not count read pairs that have their two ends mapping to  different  chromosomes
      or mapping to same chromosome but on different strands.
  donotsort:
    type: boolean?
    inputBinding:
      position: 32
      prefix: --donotsort
    doc:  |
      Do not count read pairs that have their two ends mapping to  different  chromosomes
      or mapping to same chromosome but on different strands.
  long_reads:
    type: boolean?
    inputBinding:
      position: 33
      prefix: -L
    doc:  |
      Count long reads such as Nanopore and PacBio reads. Long read counting can only run
      in one thread and  only  reads  (not  read-pairs)  can  be  counted.  There  is  no
      limitation  on  the number of 'M' operations allowed in a CIGAR string in long read
      counting.
  verbose:
    type: boolean?
    default: true
    inputBinding:
      prefix: --verbose
      position: 34
  output_file:
    type: string
    inputBinding:
      position: 35  
      prefix: -o
  inputFiles:
    type:
      type: array
      items: File
    inputBinding:
      position: 36

outputs:
  output:
    type: File
    secondaryFiles:
      - .summary
    outputBinding:
      glob: $(inputs.output_file)
  o_exon_junction_support:
    type: File?
    outputBinding:
      glob: $(inputs.output_file.concat(".jcounts"))