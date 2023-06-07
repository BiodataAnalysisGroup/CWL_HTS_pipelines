cwlVersion: v1.0
class: ExpressionTool

hints:
- class: DockerRequirement
  dockerPull: ubuntu:latest

requirements:
  InlineJavascriptRequirement: {}

inputs:
    DIRECTORY:
        type: Directory

outputs:
    raw_files: 
        type: File[]

expression: |
  ${
      var tempFile = [];

      for(var i = 0; i < inputs.DIRECTORY.listing.length; i++){
          tempFile.push(inputs.DIRECTORY.listing[i]);
      }
      
      return {
          "raw_files": tempFile
      };
  }
