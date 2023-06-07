cwlVersion: v1.0
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - '$({class: "Directory", basename: inputs.outputdir, listing: []})'

hints:
- class: DockerRequirement
  dockerPull: gawbul/docker-stringtie

baseCommand: stringtie

arguments:
- "-o"
- "$(inputs.outputdir)/$(inputs.out_gtf)"

inputs:
  mixed_reads_processing:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --mix
  long_reads_processing:
    type: boolean?
    inputBinding:
      position: 2
      prefix: -L
  expression_estimation_mode:
    type: boolean?
    inputBinding:
      position: 3
      prefix: -e
    doc:  this option directs StringTie to operate in expression estimation mode; 
          this limits the processing of read alignments to estimating the coverage of the transcripts given with the -G option 
          (hence this option requires -G).
  verbose:
    type: boolean?
    inputBinding:
      position: 4
      prefix: -v
  guide_gff:
    type: File
    inputBinding:
      prefix: -G
      position: 5
  firststrand: 
    type: boolean?
    inputBinding:
      prefix: --rf
      position: 6
    doc: Assumes a stranded library fr-firststrand.
  secondstrand: 
    type: boolean?
    inputBinding:
      prefix: --fr
      position: 7
    doc: Assumes a stranded library fr-secondstrand.
  point_features:
    type: File?
    inputBinding:
      prefix: --ptf
      position: 8
    doc: Loads a list of point-features from a text feature file <f_tab> to guide the transcriptome assembly.
  sample_label:
    type: string?
    inputBinding:
      prefix: -l
      position: 9
  min_isoform_abundance:
    type: float?
    default: 0.01
    inputBinding:
      prefix: -f
      position: 10
    doc:  Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus. 
          Lower abundance transcripts are often artifacts of incompletely spliced precursors of processed transcripts.
  cpus:
    type: int?
    inputBinding:
      prefix: -p
      position: 11
  gene_abundances:
    type: string?
    inputBinding:
      prefix: -A
      position: 12
    doc: Gene abundances will be reported (tab delimited format) in the output file with the given name.
  cov_refs:
    type: string?
    inputBinding:
      prefix: -C
      position: 13
    doc: StringTie outputs a file with the given name with all transcripts in the provided reference file that are fully covered by reads (requires -G).
  junction_filter:
    type: int?
    default: 10
    inputBinding:
      prefix: -a
      position: 14
    doc: Junctions that don't have spliced reads that align across them with at least this amount of bases on both sides are filtered out.
  junction_coverage:
    type: float?
    default: 1.0
    inputBinding: 
      prefix: -j
      position: 15
    doc:  There should be at least this many spliced reads that align across a junction (i.e. junction coverage). 
          This number can be fractional, since some reads align in more than one place. 
          A read that aligns in n places will contribute 1/n to the junction coverage.
  disable_trimming:
    type: float?
    inputBinding:
      prefix: -t
      position: 16
    doc:  This parameter disables trimming at the ends of the assembled transcripts. 
          By default StringTie adjusts the predicted transcript's start and/or stop coordinates based on sudden drops in coverage of the assembled transcript.
  min_read_coverage:
    type: float?
    default: 1.0
    inputBinding: 
      prefix: -c
      position: 17
    doc: Sets the minimum read coverage allowed for the predicted transcripts. A transcript with a lower coverage than this value is not shown in the output.
  min_read_coverage_single_exon:
    type: float?
    default: 4.75
    inputBinding: 
      prefix: -s
      position: 18
    doc: Sets the minimum read coverage allowed for single-exon transcripts.
  conservative_mode:
    type: boolean?
    inputBinding:
      prefix: --conservative
      position: 19
    doc: Assembles transcripts in a conservative mode. Same as -t -c 1.5 -f 0.05
  min_locus_gap_separation:
    type: int?
    default: 50 #(bp)
    inputBinding:
      prefix: -g
      position: 20
    doc: Minimum locus gap separation value. Reads that are mapped closer than this distance are merged together in the same processing bundle.
  ballgown_table_files:
    type: boolean?
    inputBinding: 
      prefix: -B
      position: 21
    doc:  This switch enables the output of Ballgown input table files (*.ctab) containing coverage data for the reference transcripts given with the -G option.
          With this option StringTie can be used as a direct replacement of the tablemaker program included with the Ballgown distribution.
          If the option -o is given as a full path to the output transcript file, StringTie will write the *.ctab files in the same directory as the output GTF.
  multi_location_mapped_reads:
    type: float?
    default: 0.95
    inputBinding:
      prefix: -M
      position: 22
    doc: Sets the maximum fraction of muliple-location-mapped reads that are allowed to be present at a given locus.
  ignore_read_alignments:
    type: string?
    inputBinding:
      prefix: -x
      position: 23
  turn_off_multi_mapping_correction:
    type: boolean?
    inputBinding:
      prefix: -u
      position: 24
  cram_ref:
    type: File?
    inputBinding:
      prefix: --cram-ref
      position: 25
    doc:  for CRAM input files, the reference genome sequence can be provided as a multi-FASTA file the same chromosome sequences that were used when aligning the reads.
  transcript_merge_mode:
    type: boolean?
    inputBinding: 
      prefix: --merge
      position: 26
    doc: |
      This is a special usage mode of StringTie, distinct from the assembly usage mode described above. 
      In the merge mode, StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts. 
      This mode is used in the new differential analysis pipeline to generate a global, unified set of transcripts (isoforms) across multiple RNA-Seq samples. 
      If the -G option (reference annotation) is provided, StringTie will assemble the transfrags from the input GTF files with the reference transcripts.
  out_gtf:
    type: string
  input_bam:
    type: File?
    inputBinding:
      position: 28
  outputdir:
    type: string

outputs:
  output_gtf:
    type: File
    outputBinding:
      glob: "$(inputs.outputdir)/$(inputs.out_gtf)"
  output_gene_abundances:
    type: File?
    outputBinding:
      glob: $(inputs.gene_abundances)
  output_cov_refs:
    type: File?
    outputBinding:
      glob: $(inputs.cov_refs)
  outdir:
    type: Directory?
    outputBinding:
      glob: $(inputs.outputdir)
