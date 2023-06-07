#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: kerstenbreuer/deeptools:3.1.1

baseCommand: ["multiBamSummary", "bins"]

doc: |
  Computes the read coverages for genomic regions for typically two or more BAM files. 
  The analysis can be performed for the entire genome by running the program in 'bins' mode.

inputs:
    threads: 
        type: int?
        inputBinding:
            position: 1
            prefix: -p
    labels:
        type: string[]?
        inputBinding: 
            position: 2
            prefix: -l
    smart_labels: 
        type: boolean?
        inputBinding: 
            position: 3
            prefix: --smartLabels
    binsize: 
        type: int?
        inputBinding: 
            position: 4
            prefix: -bs
    distance_between_bits:
        type: int?
        inputBinding: 
            position: 5
            prefix: -n
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
    outFileName:
        type: string?
        default: "multiBam.npz"
        inputBinding:
            position: 100
            prefix: -o
    bam:
        doc: bam file(s) as input; needs bai index file in the same directory
        type: File[]
        secondaryFiles: .bai
        inputBinding:
            position: 101
            prefix: --bamfiles

outputs:
    outNpz:
        type: File
        outputBinding:
            glob: "*.npz"
