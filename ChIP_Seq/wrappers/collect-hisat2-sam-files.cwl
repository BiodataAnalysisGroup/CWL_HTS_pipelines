cwlVersion: v1.0
class: ExpressionTool

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

expression: |
  ${
    // set variable for collecting all HISAT2 produced SAM files from both single- and paired-end fastq files
    var out = [];
    // iteration variables
    var i;
    var j;
    // iterate
    for(i = 0; i < inputs.single_files.length; i++){
      out.push(inputs.single_files[i]);
    }
    for(j = 0; j < inputs.paired_files.length; j++){
      out.push(inputs.paired_files[j]);
    }
    // return results
    return {
      "total_sam_files": out
    };
  }