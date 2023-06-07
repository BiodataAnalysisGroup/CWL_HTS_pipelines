cwlVersion: v1.0
class: ExpressionTool

hints:
- class: DockerRequirement
  dockerPull: ubuntu:latest

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: input_check
  type: string
  doc: Checks on which input files will pass to the next steps. Acceptable values are "raw_data", "fastx_trimmer" or "trim_galore".
- id: single_files
  type: File[]
- id: paired_files
  type: File[]
- id: trim_galore_single
  type: 
    type: array
    items: ["null", "File"]
- id: trim_galore_paired
  type: 
    type: array
    items: ["null", "File"]
- id: fastx_trimmer_single
  type: 
    type: array
    items: ["null", "File"]
- id: fastx_trimmer_paired
  type: 
    type: array
    items: ["null", "File"]
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
- id: single_trim
  type: File[]
- id: single_hisat2_sam
  type: string[]
- id: paired_trim_1
  type: File[]
- id: paired_trim_2
  type: File[]
- id: paired_hisat2_sam
  type: string[]

expression: |
  ${
    // set variables
    // single-end trimmed fastq
    var single_out1 = [];
    // single-end SAM filename
    var single_out2 = [];
    // paired-end trimmed fastq (forward)
    var paired_out1 = [];
    // paired-end trimmed fastq (reverse)
    var paired_out2 = [];
    // paired-end SAM filename
    var paired_out3 = [];
    // set variables to store and compare paired-end filenames with JS expressions
    var sample_ref;
    var sample_target;
    // iteration variables
    var i;
    var j;

    // fastx-trimmer
    if(inputs.input_check === "fastx_trimmer"){ // check for the presence of fastx_trimmer processed files

        // Paired-end:
        if(inputs.fastx_trimmer_paired[0] != null){ 
            for(i = 0; i < inputs.fastx_trimmer_paired.length; i++){
                if(inputs.fastx_trimmer_paired[i].basename.includes(inputs.file_split_fwd_single)){
                    // store fwd paired-end fastq file
                    paired_out1.push(inputs.fastx_trimmer_paired[i]);
                    sample_ref = inputs.fastx_trimmer_paired[i].basename.split(inputs.file_split)[0];
                    // store rev paired-end fastq file
                    for(j = 0; j < inputs.fastx_trimmer_paired.length; j++){
                        sample_target = inputs.fastx_trimmer_paired[j].basename.split(inputs.file_split)[0];
                        if(inputs.fastx_trimmer_paired[j].basename.includes(inputs.file_split_rev) && sample_ref == sample_target){
                        paired_out2.push(inputs.fastx_trimmer_paired[j]);
                        break;
                        }
                    }
                    // store SAM filename
                    paired_out3.push(sample_ref.concat(".sam"));
                }
            }
        }
        // Single-end:
        if(inputs.fastx_trimmer_single[0] != null){
            for(i = 0; i < inputs.fastx_trimmer_single.length; i++){
                // store single-end fastq file
                single_out1.push(inputs.fastx_trimmer_single[i]);
                sample_ref = inputs.fastx_trimmer_single[i].basename.split(inputs.file_split)[0];
                // store SAM filename
                single_out2.push(sample_ref.concat(".sam"));
            }        
        }

    // Trim galore
    } else if(inputs.input_check === "trim_galore"){ // check for the presence of trim_galore processed files

        // Paired-end:
        if(inputs.trim_galore_paired[0] != null){
            for(i = 0; i < inputs.trim_galore_paired.length; i++){
                if(inputs.trim_galore_paired[i].basename.includes(inputs.file_split_fwd_single)){
                    // store fwd paired-end fastq file
                    paired_out1.push(inputs.trim_galore_paired[i]);
                    sample_ref = inputs.trim_galore_paired[i].basename.split(inputs.file_split)[0];
                    // store rev paired-end fastq file
                    for(j = 0; j < inputs.trim_galore_paired.length; j++){
                        sample_target = inputs.trim_galore_paired[j].basename.split(inputs.file_split)[0];
                        if(inputs.trim_galore_paired[j].basename.includes(inputs.file_split_rev) && sample_ref == sample_target){
                        paired_out2.push(inputs.trim_galore_paired[j]);
                        break;
                        }
                    }
                    // store SAM filename
                    paired_out3.push(sample_ref.concat(".sam"));
                }
            }
        }
        // Single-end:
        if(inputs.trim_galore_single[0] != null){
            for(i = 0; i < inputs.trim_galore_single.length; i++){
                // store single-end fastq file
                single_out1.push(inputs.trim_galore_single[i]);
                sample_ref = inputs.trim_galore_single[i].basename.split(inputs.file_split)[0];
                // store SAM filename
                single_out2.push(sample_ref.concat(".sam"));
            }
        }

      // Raw data
    } else if(inputs.input_check === "raw_data"){ // finally check if trimming is set to false and raw files should be used

      // Paired-end:
      for(i = 0; i < inputs.paired_files.length; i++){
        if(inputs.paired_files[i].basename.includes(inputs.file_split_fwd_single)){
            // store fwd paired-end fastq file
            paired_out1.push(inputs.paired_files[i]);
            sample_ref = inputs.paired_files[i].basename.split(inputs.file_split)[0];
            // store rev paired-end fastq file
            for(j = 0; j < inputs.paired_files.length; j++){
                sample_target = inputs.paired_files[j].basename.split(inputs.file_split)[0];
                if(inputs.paired_files[j].basename.includes(inputs.file_split_rev) && sample_ref == sample_target){
                paired_out2.push(inputs.paired_files[j]);
                break;
                }
            }
            // store SAM filename
            paired_out3.push(sample_ref.concat(".sam"));
        }
      }  

      // Single-end:
      for(i = 0; i < inputs.single_files.length; i++){
        // store single-end fastq file
        single_out1.push(inputs.single_files[i]);
        sample_ref = inputs.single_files[i].basename.split(inputs.file_split)[0];
        // store SAM filename
        single_out2.push(sample_ref.concat(".sam"));
      }  

    }
    // return results
    return {
        "single_trim": single_out1, 
        "single_hisat2_sam": single_out2, 
        "paired_trim_1": paired_out1, 
        "paired_trim_2": paired_out2,
        "paired_hisat2_sam": paired_out3
    };
  }
