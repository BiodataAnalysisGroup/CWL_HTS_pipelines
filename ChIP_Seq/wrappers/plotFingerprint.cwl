#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: kerstenbreuer/deeptools:3.1.1

baseCommand: ["plotFingerprint"]
  
inputs:
    bam:
        doc: bam file(s) as input; needs bai index file in the same directory
        type: File[]
        secondaryFiles: .bai
        inputBinding:
            position: 100
            prefix: --bamfiles
    threads: 
        type: int?
        default: 8
        inputBinding:
            position: 1
            prefix: -p
    outFileName:
        type: string?
        default: "fingerprint_plot.png"
        inputBinding:
            position: 2
            prefix: -plot
    plotFileFormat: 
        type: string?
        inputBinding: 
            position: 3
            prefix: --plotFileFormat
    blackListFileName: 
        type: File? 
        inputBinding: 
            position: 4
            prefix: -bl
        doc: |
            A BED or GTF file containing regions that should be excluded from all analyses. 
            Currently this works by rejecting genomic chunks that happen to overlap an entry. 
            Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, 
            then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.

outputs:
    outImage:
        type: File
        outputBinding:
            glob: $(inputs.outFileName)
