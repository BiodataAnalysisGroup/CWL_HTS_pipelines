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
- id: input_paired
  type: File[]
- id: trimming_single
  type: File[]
- id: trimming_paired
  type: File[]

outputs:
- id: single_files
  type: File[]
- id: paired_files
  type: File[]

expression: |
  ${
    var single;
    var paired;
    //
    if(inputs.trimming_check){
        //
        single = inputs.trimming_single;
        paired = inputs.trimming_paired;
    } else {
        //
        single = inputs.input_single;
        paired = inputs.input_paired;
    }
    //
    return {
        "single_files": single,
        "paired_files": paired
    };
  }