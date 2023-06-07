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
- id: qc_check
  type: boolean
- id: trimming_check
  type: boolean

outputs:
- id: single_files
  type: File[]
- id: paired_files
  type: File[]
- id: trim_galore_for_single
  type: string
- id: trim_galore_for_paired
  type: string
- id: fastqc_for_raw
  type: string
- id: fastqc_for_single
  type: string
- id: fastqc_for_paired
  type: string
- id: cp_command_raw
  type: string
- id: cp_command_single
  type: string
- id: cp_command_paired
  type: string

expression: |
  ${
    var single = [];
    var paired = [];
    var flag_single;
    var base_ref;
    var base_target;

    var single_trim_command = "/usr/local/bin/trim_galore";
    var paired_trim_command = "/usr/local/bin/trim_galore";

    var raw_fastqc_command = "/opt/FastQC/fastqc";
    var single_fastqc_command = "/opt/FastQC/fastqc";
    var paired_fastqc_command = "/opt/FastQC/fastqc";

    var command_raw = "cp"
    var command_single = "cp"
    var command_paired = "cp"

    for(var i = 0; i < inputs.input_files.length; i++){

      if(inputs.input_files[i].basename.includes(inputs.file_split_fwd_single)){
        base_ref = inputs.input_files[i].basename.split(inputs.file_split)[0];

        flag_single = 1;

        for(var j = 0; j < inputs.input_files.length; j++){
          base_target = inputs.input_files[j].basename.split(inputs.file_split)[0];

          if(base_ref == base_target && inputs.input_files[j].basename.includes(inputs.file_split_rev)){
            flag_single = 0;
            paired.push(inputs.input_files[i]);
            paired.push(inputs.input_files[j]);
            break;
          }
        }

        if(flag_single == 1){
          single.push(inputs.input_files[i]);
        }
      }
      
    }

    if(single.length == 0){
      single_trim_command = "echo";
      single_fastqc_command = "echo";
      //
      command_single = "echo";
    }

    if(paired.length == 0){
      paired_trim_command = "echo";
      paired_fastqc_command = "echo";
      //
      command_paired = "echo";
    }

    // QC and trimming check
    if(!inputs.qc_check){
        //
        raw_fastqc_command = "echo";
        single_fastqc_command = "echo";
        paired_fastqc_command = "echo";
        command_raw = "echo";
        command_single = "echo";
        command_paired = "echo";
    }
    //
    if(!inputs.trimming_check){
        //
        single_trim_command = "echo";
        paired_trim_command = "echo";
    }

    return {
      "single_files": single, 
      "paired_files": paired,
      "trim_galore_for_single": single_trim_command,
      "trim_galore_for_paired": paired_trim_command,
      "fastqc_for_raw": raw_fastqc_command,
      "fastqc_for_single": single_fastqc_command,
      "fastqc_for_paired": paired_fastqc_command,
      "cp_command_raw": command_raw,
      "cp_command_single": command_single,
      "cp_command_paired": command_paired
    };
  }