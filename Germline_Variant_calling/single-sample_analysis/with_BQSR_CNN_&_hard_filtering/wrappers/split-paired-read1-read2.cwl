cwlVersion: v1.0
class: ExpressionTool

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: paired_files
  type: File[]
- id: file_split
  type: string
  default: "_R"
- id: file_split_fwd_single
  type: string
  default: "R1"
- id: file_split_rev
  type: string
  default: "R2"

outputs:
- id: reads_1
  type: File[]
- id: reads_2
  type: File[]

expression: |
  ${
    // paired-end fwd and rev file variables
    var read1 = [];
    var read2 = [];
    // variables for filename comparisons with JS expressions
    var base_ref;
    var base_target;
    // iterate
    for(var i = 0; i < inputs.paired_files.length; i++){
      if(inputs.paired_files[i].basename.includes(inputs.file_split_fwd_single)){
        // detect fwd paired-end filename
        base_ref = inputs.paired_files[i].basename.split(inputs.file_split)[0];
        // iterate
        for(var j = 0; j < inputs.paired_files.length; j++){
          base_target = inputs.paired_files[j].basename.split(inputs.file_split)[0];
          // detect rev paired-end filename and save
          if(base_ref == base_target && inputs.paired_files[j].basename.includes(inputs.file_split_rev)){
            read1.push(inputs.paired_files[i]);
            read2.push(inputs.paired_files[j]);
            break;
          }
        }
        
      }
    }
    // return results
    return {
      "reads_1": read1, 
      "reads_2": read2
    };
  }