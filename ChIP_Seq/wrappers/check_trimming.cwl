cwlVersion: v1.0
class: ExpressionTool

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: file_split
  type: string
  default: "_R"
- id: file_split_fwd_single
  type: string
  default: "R1"
- id: file_split_rev
  type: string
  default: "R2"
- id: trimming_check
  type: boolean
- id: input_single
  type: File[]
- id: input_paired_fwd
  type: File[]
- id: input_paired_rev
  type: File[]
- id: trimming_single
  type: 
  - type: array
    items: ["null", File]
- id: trimming_paired_fwd
  type: 
  - type: array
    items: ["null", File]
- id: trimming_paired_rev
  type: 
  - type: array
    items: ["null", File]

outputs:
- id: single_files
  type: File[]
- id: paired_files_fwd
  type: File[]
- id: paired_files_rev
  type: File[]

expression: |
  ${
    // set variables
    var single;
    var paired_fwd;
    var paired_rev;
    // check user option for trimming
    if(inputs.trimming_check){
        single = inputs.trimming_single;
        paired_fwd = inputs.trimming_paired_fwd;
        paired_rev = inputs.trimming_paired_rev;
    } else {
        single = inputs.input_single;
        paired_fwd = inputs.input_paired_fwd;
        paired_rev = inputs.input_paired_rev;
    }
    // return results
    return {
        "single_files": single,
        "paired_files_fwd": paired_fwd,
        "paired_files_rev": paired_rev
    };
  }