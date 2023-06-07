cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: broadinstitute/gatk:4.3.0.0

baseCommand: [gatk, CombineGVCFs]

arguments:
- shellQuote: false
  position: 1
  valueFrom: |
   ${
    var args = [];
    for(var i = 0; i < inputs.gvcf_files.length; i++){
        var arg = inputs.gvcf_files[i].path;
        args.push("-V " + arg);
    }
    var res = args.join(" ");
    return res;
   }

inputs:
    reference: 
        type: File
        secondaryFiles:
        - .amb
        - .ann
        - .bwt
        - .pac
        - .sa
        - .fai
        - ^.dict
        inputBinding:
            position: 2
            prefix: -R
    gvcf_files:
        type: File[]
        secondaryFiles: 
        - .tbi
    intervals:
        type: File?
        inputBinding: 
            position: 3
            prefix: -L
    exclude_intervals:
        type: File?
        inputBinding:
            position: 4
            prefix: -XL
    output_name:
        type: string
        default: "cohort.g.vcf.gz"
        inputBinding: 
            position: 5
            prefix: -O

outputs: 
    output:
        type: File
        secondaryFiles:
        - .tbi
        outputBinding:
            glob: $(inputs.output_name)
