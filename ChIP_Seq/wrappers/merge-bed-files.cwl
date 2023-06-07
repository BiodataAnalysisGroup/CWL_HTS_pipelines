cwlVersion: v1.0
class: ExpressionTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: ubuntu:latest

inputs:
- id: header_file
  type: File
- id: bed_files
  type: File[]

outputs:
- id: total_bed_files
  type: File[]

expression: |
  ${
    var out = [];
    var i;
    
    out.push(inputs.header_file)

    for(i = 0; i < inputs.bed_files.length; i++){
      out.push(inputs.bed_files[i]);
    }

    return {
      "total_bed_files": out
    };
  }
