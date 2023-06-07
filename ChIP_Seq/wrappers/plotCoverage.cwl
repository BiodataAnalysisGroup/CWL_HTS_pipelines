#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: kerstenbreuer/deeptools:3.1.1

doc: |
  Create a plot showing global coverage for reads provided in BAM.

baseCommand: ["plotCoverage"]
  
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
    skipZeros:
        type: boolean?
        inputBinding:
            position: 2
            prefix: --skipZeros
    verbose:
        type: boolean?
        inputBinding:
            position: 3
            prefix: --verbose
    plotFileFormat: 
        type: string?
        inputBinding: 
            position: 4
            prefix: --plotFileFormat
    outFileName:
        type: string
        default: "coverage_plot.png"
        inputBinding:
            position: 5
            prefix: -o
    region: 
        type: string?
        inputBinding: 
            position: 6
            prefix: -r
    blackListFileName: 
        type: File?
        inputBinding: 
            position: 7
            prefix: -bl
    bed_files: 
        type: File[]?
        inputBinding: 
            position: 8
            prefix: --BED

outputs:
  outImage:
        type: File
        outputBinding:
            glob: $(inputs.outFileName)
