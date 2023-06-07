cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: kerstenbreuer/deeptools:3.1.1

doc: |
  This tool calculates scores per genome regions and prepares an intermediate file that can be used with plotHeatmap and plotProfiles.

baseCommand: ["computeMatrix"]

arguments: 
    - "reference-point"

inputs:
    bw:
        doc: bigwig files as input
        type: File[]
        inputBinding:
            position: 101
            prefix: -S
    regions: 
        doc: File name or names, in BED or GTF format, containing the regions to plot. 
        type: File
        inputBinding: 
            position: 100
            prefix: -R
    upstream:
        doc: Distance upstream of the start site of the regions defined in the region file. If the regions are genes, this would be the distance upstream of the transcription start site.
        type: int?
        inputBinding:
            position: 1
            prefix: -b
    downstream: 
        doc: Distance downstream of the end site of the given regions. If the regions are genes, this would be the distance downstream of the transcription end site.
        type: int?
        inputBinding:
            position: 2
            prefix: -a
    skipZeros:
        type: boolean?
        default: true
        inputBinding:
            position: 3
            prefix: --skipZeros
    outputFile: 
        type: string?
        default: "matrix.out"
        inputBinding:
            position: 4
            prefix: -o
    outFileSortedRegions:
        type: string?
        default: "regions.out"
        inputBinding:
            position: 5
            prefix: --outFileSortedRegions
    threads: 
        type: int?
        default: 8
        inputBinding:
            position: 6
            prefix: -p
    blackListFileName: 
        type: File? 
        inputBinding: 
            position: 9
            prefix: -bl
        doc: |
            A BED or GTF file containing regions that should be excluded from all analyses. 
            Currently this works by rejecting genomic chunks that happen to overlap an entry. 
            Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, 
            then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant.

outputs:
    outMatrix:
        type: File
        outputBinding:
            glob: $(inputs.outputFile)
    outRegions:
        type: File
        outputBinding:
            glob: $(inputs.outFileSortedRegions)
