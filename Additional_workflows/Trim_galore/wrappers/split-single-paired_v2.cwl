cwlVersion: v1.0
class: ExpressionTool

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: input_files
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
- id: single_files
  type: File[]
- id: paired_files_fwd
  type: File[]
- id: paired_files_rev
  type: File[]
# strings for trimming commands
- id: trim_galore_command_single
  type: string
- id: trim_galore_command_paired
  type: string  

expression: |
  ${
    // define variables
    var single = [];
    var paired_fwd = [];
    var paired_rev = [];
    // trimmomatic commands
    var single_trim_command = "trim_galore";
    var paired_trim_command = "trim_galore";
    // iterate through files
    for(var i = 0; i < inputs.input_files.length; i++){
      // check for single- and paired-end forward files 
      if(inputs.input_files[i].basename.includes(inputs.file_split_fwd_single)){
        base_ref = inputs.input_files[i].basename.split(inputs.file_split)[0];
        flag_single = 1;
        for(var j = 0; j < inputs.input_files.length; j++){
          base_target = inputs.input_files[j].basename.split(inputs.file_split)[0];
          // check for paired-end reverse files
          if(base_ref == base_target && inputs.input_files[j].basename.includes(inputs.file_split_rev)){
            flag_single = 0;
            paired_fwd.push(inputs.input_files[i]);
            paired_rev.push(inputs.input_files[j]);
            break;
          }
        }
        // if single-end files exist
        if(flag_single == 1){
          single.push(inputs.input_files[i]);
        }
      }
    }
    // if single-end files do not exist
    if(single.length == 0){
      single_trim_command = "echo";
    }
    // if paired-end files do not exist
    if(paired_fwd.length == 0){
      paired_trim_command = "echo";
    }
    // return results
    return {
      // total files
      "single_files": single, 
      "paired_files_fwd": paired_fwd,
      "paired_files_rev": paired_rev,
      // trimming commands
      "trim_galore_command_single": single_trim_command,
      "trim_galore_command_paired": paired_trim_command
    };
  }