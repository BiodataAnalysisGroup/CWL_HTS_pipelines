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
    // set variables
    // single-end fastq files
    var single = [];
    // paired-end fastq files
    var paired = [];
    // variables for comparing filenames and separating them based on JS expressions
    var flag_single;
    var base_ref;
    var base_target;
    // default command for trimming (stays the same or changes based on user input)
    // trim galore
    var single_trim_command = "/usr/local/bin/trim_galore";
    var paired_trim_command = "/usr/local/bin/trim_galore";

    // default command for quality control (stays the same or changes based on user input)
    // fastqc
    var raw_fastqc_command = "/opt/FastQC/fastqc";
    var single_fastqc_command = "/opt/FastQC/fastqc";
    var paired_fastqc_command = "/opt/FastQC/fastqc";
    // cp
    var command_raw = "cp"
    var command_single = "cp"
    var command_paired = "cp"

    // iterate through list of input fastq files
    for(var i = 0; i < inputs.input_files.length; i++){
      // check for single- or paired end files based on file_split user input
      if(inputs.input_files[i].basename.includes(inputs.file_split_fwd_single)){
        base_ref = inputs.input_files[i].basename.split(inputs.file_split)[0];
        // single-end or paired-end fwd fastq file detected
        flag_single = 1;
        // check for the presence of rev paired-end fastq file
        for(var j = 0; j < inputs.input_files.length; j++){
          base_target = inputs.input_files[j].basename.split(inputs.file_split)[0];
          // store paired-end fastq files and set flag_single = 0
          if(base_ref == base_target && inputs.input_files[j].basename.includes(inputs.file_split_rev)){
            flag_single = 0;
            paired.push(inputs.input_files[i]);
            paired.push(inputs.input_files[j]);
            break;
          }
        }
        // if flag_single = 1 then no rev paired-end fastq file was found
        // store single-end fastq file
        if(flag_single == 1){
          single.push(inputs.input_files[i]);
        }
      }
      
    }
    // define single- and paired-end QC and trimming commands based on the detected files
    // single-end pre-processing commands
    if(single.length == 0){
      single_trim_command = "echo";
      single_fastqc_command = "echo";
      command_single = "echo";
    }
    // paired-end pre-processing commands
    if(paired.length == 0){
      paired_trim_command = "echo";
      paired_fastqc_command = "echo";
      command_paired = "echo";
    }

    // check user options for QC and trimming and define the respective commands
    // QC
    if(!inputs.qc_check){
        raw_fastqc_command = "echo";
        single_fastqc_command = "echo";
        paired_fastqc_command = "echo";
        command_raw = "echo";
        command_single = "echo";
        command_paired = "echo";
    }
    // trimming
    if(!inputs.trimming_check){
        single_trim_command = "echo";
        paired_trim_command = "echo";
    }
    // return results
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