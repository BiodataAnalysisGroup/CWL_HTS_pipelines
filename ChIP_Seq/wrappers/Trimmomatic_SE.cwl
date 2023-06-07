cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: "staphb/trimmomatic:latest"

arguments: 
- position: 1
  shellQuote: false
  valueFrom: $(inputs.command)
- position: 2
  shellQuote: false
  valueFrom: "SE"
- position: 8
  valueFrom: $( inputs.input_fastq.basename.split(inputs.file_split)[0].concat(inputs.file_split_fwd_single, "_trimmed.fastq.gz") )

inputs:
    command:
        type: string
        default: "trimmomatic"
    file_split:
        type: string
        default: "_R"
    file_split_fwd_single:
        type: string
        default: "R1"
    trimm_se_threads:
        type: int?
        default: 8
        inputBinding:
            position: 3
            prefix: -threads
    phred33:
        type: boolean? 
        inputBinding: 
            position: 4
            prefix: -phred33
    phred64:
        type: boolean? 
        inputBinding: 
            position: 5
            prefix: -phred64
    trimLogStr:
        type: string?
        inputBinding:
            position: 6
            prefix: -trimlog
    input_fastq:
        type: File
        inputBinding:
            position: 7
    # trimming parameters
    illuminaClip: 
        type: string?
        inputBinding:
            position: 9
            prefix: "ILLUMINACLIP:"
            separate: false
            shellQuote: false
    slidingWindow:
        type: string?
        default: "4:18"
        inputBinding: 
            position: 10
            prefix: "SLIDINGWINDOW:"
            separate: false
            shellQuote: false
    leading:
        type: int?
        default: 28
        inputBinding: 
            position: 11
            prefix: "LEADING:" 
            separate: false
            shellQuote: false
    trailing:
        type: int?
        default: 28
        inputBinding: 
            position: 12
            prefix: "TRAILING:" 
            separate: false
            shellQuote: false
    minlen:
        type: int?
        default: 36
        inputBinding: 
            position: 13
            prefix: "MINLEN:" 
            separate: false
            shellQuote: false
    maxinfo: 
        type: string?
        inputBinding:
            position: 14
            prefix: "MAXINFO:"
            separate: false
            shellQuote: false
    crop:
        type: int?
        inputBinding: 
            position: 15
            prefix: "CROP:"
            separate: false
            shellQuote: false
    headcrop: 
        type: int?
        inputBinding: 
            position: 16
            prefix: "HEADCROP:"
            separate: false
            shellQuote: false

stderr: $( inputs.input_fastq.basename.split(inputs.file_split)[0].concat('_trimmed.log') )

outputs:
    stderr_log:
        type: stderr
    outFastq:
        type: File?
        outputBinding:
            glob: $( inputs.input_fastq.basename.split(inputs.file_split)[0].concat(inputs.file_split_fwd_single, "_trimmed.fastq.gz") ) 
            outputEval: 
             ${
                if(self.length === 0){
                    return {};
                } else {
                    return self;
                }
             }
    trimLog_out:
        type: File?
        outputBinding:
            glob: $(inputs.trimLogStr)
