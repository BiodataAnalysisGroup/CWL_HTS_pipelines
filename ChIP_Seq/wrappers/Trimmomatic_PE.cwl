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
  valueFrom: "PE"
- position: 9
  valueFrom: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat(inputs.file_split_fwd_single, "_trimmed.fastq.gz") )
- position: 10
  valueFrom: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat(inputs.file_split_fwd_single, "_trimmed_unpaired.fastq.gz") )
- position: 11
  valueFrom: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat(inputs.file_split_rev, "_trimmed.fastq.gz") )
- position: 12
  valueFrom: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat(inputs.file_split_rev, "_trimmed_unpaired.fastq.gz") )

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
    file_split_rev:
        type: string
        default: "R2"
    trimm_pe_threads:
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
    # Input paired-end fastq files
    input_fastq_fwd:
        type: File
        inputBinding:
            position: 7
    input_fastq_rev:
        type: File
        inputBinding:
            position: 8
    # trimming parameters
    illuminaClip: 
        type: string?
        inputBinding:
            position: 13
            prefix: "ILLUMINACLIP:"
            separate: false
            shellQuote: false
    slidingWindow:
        type: string
        default: "4:18"
        inputBinding: 
            position: 14
            prefix: "SLIDINGWINDOW:"
            separate: false
            shellQuote: false
    leading:
        type: int
        default: 28
        inputBinding: 
            position: 15
            prefix: "LEADING:" 
            separate: false
            shellQuote: false
    trailing:
        type: int
        default: 28
        inputBinding: 
            position: 16
            prefix: "TRAILING:" 
            separate: false
            shellQuote: false
    minlen:
        type: int
        default: 36
        inputBinding: 
            position: 17
            prefix: "MINLEN:" 
            separate: false
            shellQuote: false
    maxinfo: 
        type: string?
        inputBinding:
            position: 18
            prefix: "MAXINFO:"
            separate: false
            shellQuote: false
    crop:
        type: int?
        inputBinding: 
            position: 19
            prefix: "CROP:"
            separate: false
            shellQuote: false
    headcrop: 
        type: int?
        inputBinding: 
            position: 20
            prefix: "HEADCROP:"
            separate: false
            shellQuote: false

stderr: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat('_trimmed.log') )

outputs:
    # stderr log
    stderr_log:
        type: stderr
    # output FASTQ files
    outFastq_fwd_paired:
        type: File?
        outputBinding:
            glob: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat(inputs.file_split_fwd_single, "_trimmed.fastq.gz") ) 
            outputEval: 
             ${
                if(self.length === 0){
                    return {};
                } else {
                    return self;
                }
             }
    outFastq_fwd_unpaired:
        type: File?
        outputBinding:
            glob: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat(inputs.file_split_fwd_single, "_trimmed_unpaired.fastq.gz") ) 
            outputEval: 
             ${
                if(self.length === 0){
                    return {};
                } else {
                    return self;
                }
             }
    outFastq_rev_paired:
        type: File?
        outputBinding:
            glob: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat(inputs.file_split_rev, "_trimmed.fastq.gz") ) 
            outputEval: 
             ${
                if(self.length === 0){
                    return {};
                } else {
                    return self;
                }
             }
    outFastq_rev_unpaired:
        type: File?
        outputBinding:
            glob: $( inputs.input_fastq_fwd.basename.split(inputs.file_split)[0].concat(inputs.file_split_rev, "_trimmed_unpaired.fastq.gz") ) 
            outputEval: 
             ${
                if(self.length === 0){
                    return {};
                } else {
                    return self;
                }
             }
    # trimmomatic trimming log
    trimLog_out:
        type: File?
        outputBinding:
            glob: $(inputs.trimLogStr)
