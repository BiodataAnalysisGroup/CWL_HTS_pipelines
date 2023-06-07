#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: DockerRequirement
    dockerPull: kerstenbreuer/deeptools:3.1.1

doc: |
  Generates coverage tracks in bedgraph or bigiwig format from BAM file.
  Normalization by spike-in reads is supported.

baseCommand: ["bamCoverage"]

inputs:
  bam:
    doc: bam file as input; needs bai index file in the same directory
    type: File
    secondaryFiles: .bai
    inputBinding:
        position: 100
        prefix: --bam
  threads: 
        type: int?
        default: 8
        inputBinding:
            position: 2
            prefix: -p
  effective_genome_size:
    doc: |
      the effectively mappable genome size, 
      see: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
    type: long?
    inputBinding:
        position: 3
        prefix: --effectiveGenomeSize
  bin_size:
    type: int?
    inputBinding:
      prefix: --binSize
      position: 4
  ignoreForNormalization:
    type: string[]?
    inputBinding:
      prefix: --ignoreForNormalization
      position: 5
  normalizeUsing:
    type: string
  extendReads:
    doc: |
      This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception. 
      NOTE: This feature is generally NOT recommended for spliced-read data, such as RNA-seq, as it would extend reads over skipped regions. 
      Single-end: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended. 
      Paired-end: Reads with mates are always extended to match the fragment size defined by the two read mates. 
      Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. 
      The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads).
    type: int?
    inputBinding:
      prefix: --extendReads
      position: 6
  spike_in_count:
    doc: |
      Number of reads aligned to the spike in reference, optional.
      If specified, coverage will be multiplied by 1/spike_in_count and
      normalizeUsing will be ignored.
    type: long?
    inputBinding:
      position: 7
      prefix: --scaleFactor
      valueFrom: |
        ${ 
          if( self == null ){
            return null
          }
          else{
            return (1.0 / parseFloat(self)) 
          }
        }
  outFileFormat:
    doc: bigwig or bedgraph
    type: string
    default: bigwig
    inputBinding:
      position: 8
      prefix: --outFileFormat
  blackListFileName: 
    type: File? 
    inputBinding: 
      position: 9
      prefix: -bl
    doc: |
     A BED or GTF file containing regions that should be excluded from all analyses. 
     Currently this works by rejecting genomic chunks that happen to overlap an entry. 
     Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, 
     then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.

arguments:  
  - valueFrom: $(inputs.bam.nameroot).bw
    prefix: --outFileName
    position: 9
  - valueFrom: |
      ${ 
        if( inputs.spike_in_count == null ){
          return inputs.normalizeUsing
        }
        else{
          return null 
        }
      }
    prefix: --normalizeUsing
    position: 10

outputs:
  bigwig:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot).bw
