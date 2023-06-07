cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement


baseCommand: [printf]

arguments: 
- position: 1
  valueFrom: |
   ${
    var fileNames = ["chromosome", "start", "end"];
    var i;
    // iterate the array of files to extract filenames
    for(i = 0; i < inputs.input_files.length; i++){
        // add file
        fileNames.push(inputs.input_files[i].basename);
    }
    return fileNames.join("\t");
   }

inputs:
    input_files:
        type: File[]
    output_name:
        type: string

outputs:
    output:
        type: stdout

stdout: $(inputs.output_name)
