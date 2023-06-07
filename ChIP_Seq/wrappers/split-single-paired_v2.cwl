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
- id: paired_files_fwd
  type: File[]
- id: paired_files_rev
  type: File[]
# strings for fastq file basenames
- id: single_files_basenames
  type: string[]
- id: paired_files_basenames
  type: string[]
# strings for fastqc commands
- id: fastqc_for_raw
  type: string
- id: fastqc_for_single
  type: string
- id: fastqc_for_paired
  type: string
# strings for cp commands
- id: cp_command_raw
  type: string
- id: cp_command_single
  type: string
- id: cp_command_paired
  type: string
# strings for trimming commands
- id: trimmomatic_command_single
  type: string
- id: trimmomatic_command_paired
  type: string  
# string for single- and paired-end basenames
- id: single_files_sam
  type: string[]
- id: paired_files_sam
  type: string[]

expression: |
  ${
    // define variables
    var single = [];
    var paired_fwd = [];
    var paired_rev = [];
    var flag_single;
    var base_ref;
    var base_target;
    // total fastq files
    var single_names = [];
    var paired_names = [];
    // SAM files
    var single_sam = [];
    var paired_sam = [];
    // trimmomatic commands
    var single_trim_command = "trimmomatic";
    var paired_trim_command = "trimmomatic";
    // fastqc commands
    var raw_fastqc_command = "/opt/FastQC/fastqc";
    var single_fastqc_command = "/opt/FastQC/fastqc";
    var paired_fastqc_command = "/opt/FastQC/fastqc";
    // cp commands
    var command_raw = "cp";
    var command_single = "cp";
    var command_paired = "cp";
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
            paired_names.push( inputs.input_files[i].basename.split(inputs.file_split)[0] );
            paired_sam.push( inputs.input_files[i].basename.split(inputs.file_split)[0].concat(".sam") );
            break;
          }
        }
        // if single-end files exist
        if(flag_single == 1){
          single.push(inputs.input_files[i]);
          single_names.push( inputs.input_files[i].basename.split(inputs.file_split)[0] );
          single_sam.push( inputs.input_files[i].basename.split(inputs.file_split)[0].concat(".sam") );
        }
      }
    }
    // if single-end files do not exist
    if(single.length == 0){
      single_trim_command = "echo";
      single_fastqc_command = "echo";
      command_single = "echo";
    }
    // if paired-end files do not exist
    if(paired_fwd.length == 0){
      paired_trim_command = "echo";
      paired_fastqc_command = "echo";
      command_paired = "echo";
    }
    // QC option check
    if(!inputs.qc_check){
        raw_fastqc_command = "echo";
        single_fastqc_command = "echo";
        paired_fastqc_command = "echo";
        command_raw = "echo";
        command_single = "echo";
        command_paired = "echo";
    }
    // trimming option check
    if(!inputs.trimming_check){
        single_trim_command = "echo";
        paired_trim_command = "echo";
    }
    // return results
    return {
      // total files
      "single_files": single, 
      "paired_files_fwd": paired_fwd,
      "paired_files_rev": paired_rev,
      // basenames
      "single_files_basenames": single_names,
      "paired_files_basenames": paired_names,
      // SAM files
      "single_files_sam": single_sam,
      "paired_files_sam": paired_sam,
      // FASTQC commands
      "fastqc_for_raw": raw_fastqc_command,
      "fastqc_for_single": single_fastqc_command,
      "fastqc_for_paired": paired_fastqc_command,
      // cp commands
      "cp_command_raw": command_raw,
      "cp_command_single": command_single,
      "cp_command_paired": command_paired,
      // trimmomatic commands
      "trimmomatic_command_single": single_trim_command,
      "trimmomatic_command_paired": paired_trim_command
    };
  }