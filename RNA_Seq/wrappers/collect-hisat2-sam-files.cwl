cwlVersion: v1.0
class: ExpressionTool

hints:
- class: DockerRequirement
  dockerPull: ubuntu:latest

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: single_files
  type: File[]
- id: paired_files
  type: File[]

outputs:
- id: total_sam_files
  type: File[]
- id: names_samtools_view_unsorted
  type: string[]
- id: names_samtools_sort_sorted
  type: string[]
- id: names_stringtie_gtf
  type: string[]
- id: names_htseq_counts
  type: string[]
- id: names_basenames
  type: string[]

expression: |
  ${
    // set variable for collecting all HISAT2 produced SAM files from both single- and paired-end fastq files
    var out = [];
    // return results
    return {
      "total_sam_files": out
    };
  }