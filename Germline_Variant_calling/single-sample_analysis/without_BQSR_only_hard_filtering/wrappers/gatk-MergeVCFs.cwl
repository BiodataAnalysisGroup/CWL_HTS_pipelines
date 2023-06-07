cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, MergeVcfs]

arguments:
- shellQuote: false
  prefix: -O
  valueFrom: $(inputs.OUTPUT + ".vcf.gz")
- shellQuote: false
  valueFrom: |
   ${
    var args = [];
    for(var i = 0; i < inputs.INPUT.length; i++){
        var arg = inputs.INPUT[i].path;
        args.push("-I " + arg);
    }
    var res = args.join(" ");
    return res;
   }

inputs:
  INPUT: 
    type: File[]
  OUTPUT: 
    type: string

outputs:
  output:
    type: File
    secondaryFiles:
    - .tbi
    outputBinding: 
      glob: $(inputs.OUTPUT + ".vcf.gz")
