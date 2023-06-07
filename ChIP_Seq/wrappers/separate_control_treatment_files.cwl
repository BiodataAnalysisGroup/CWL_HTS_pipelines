cwlVersion: v1.0
class: ExpressionTool

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: treatment_samples
  type: string[]
  doc: |
   A list of treatment sample names for MACS2 and differential binding analysis.
- id: control_samples
  type: string[]?
  doc: |
   A list of control sample names for MACS2 and differential binding analysis. The same control sample name is repeated for its different treatment samples. 
- id: aligned_files
  type: File[]
  secondaryFiles: .bai

outputs:
- id: treatment_files
  type: File[]
  secondaryFiles: .bai
- id: control_files
  type: File[]?
  secondaryFiles: .bai

expression: |
 ${
    // var for trearment samples
    var trm = [];
    // var for control samples
    var ctl = [];
    // iteration var
    var i;
    var j;
    // iterate through treatment and control sample names (arrays of equal length)
    for(i = 0; i < inputs.treatment_samples.length; i++){
        // iterate through aligned files (BAM)
        for(j = 0; j < inputs.aligned_files.length; j++){
            // add treatment file
            if(inputs.aligned_files[j].basename.includes(inputs.treatment_samples[i])){
                trm.push(inputs.aligned_files[j]);
            }
            // check for control sample names
            if(inputs.control_samples.length !== 0){
                // add control file
                if(inputs.aligned_files[j].basename.includes(inputs.control_samples[i])){
                    ctl.push(inputs.aligned_files[j]);
                }
            } else {
              ctl.push(null);
            }
        }
    }
    // return arrays of files
    return {
        "treatment_files": trm,
        "control_files": ctl
    };
 }
