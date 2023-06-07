cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: DockerRequirement
  dockerPull: broadinstitute/gatk:4.3.0.0

arguments:
    - gatk
    - GatherBQSRReports
    - prefix: "-O"
      valueFrom: $(inputs.output_file + ".bqsr.table")
    - shellQuote: false
      valueFrom: | 
        ${
            var args = [];
            for(var i = 0; i < inputs.bqsr_tables.length; i++){
                var arg = inputs.bqsr_tables[i].path;
                args.push("-I " + arg);
            }
            var res = args.join(" ");
            return res;
        }
    
inputs:
    input: 
        type: File?
        inputBinding: 
            position: 1
            prefix: -I
    arguments_file:
        type: File?
        inputBinding: 
            position: 2
            prefix: --arguments_file
    bqsr_tables:
        type: File[]?
    output_file: 
        type: string

outputs:
    output: 
        type: File
        outputBinding: 
            glob: $(inputs.output_file + ".bqsr.table")
