cwlVersion: v1.0
class: CommandLineTool

requirements:
    InlineJavascriptRequirement: {}

hints:
- class: DockerRequirement
  dockerPull: ewels/multiqc:v1.13

# Main input(s) is the directory containing fastqc .zip reports

inputs:
    force: 
        type: boolean?
        default: false
        inputBinding:
            position: 1
            prefix: -f
        doc: Overwrite any existing reports
    input_dir:
        type: Directory?
        inputBinding: 
            position: 2
            prefix: -d
        doc: Prepend directory to sample names
    fullnames: 
        type: boolean?
        default: false
        inputBinding:
            position: 3
            prefix: -s
        doc: Do not clean the sample names (leave as full file name)
    title: 
        type: string?
        inputBinding:
            position: 4
            prefix: -i
        doc: Report title. Printed as page header, used for filename if not otherwise specified.
    filename: 
        type: string?
        default: "multiqc"
        inputBinding:
            position: 5
            prefix: -n
        doc: Report filename. Use 'stdout' to print to standard out.
    verbose: 
        type: boolean?
        default: true
        inputBinding:
            position: 6
            prefix: -v
        doc: Increase output verbosity.
    filelist:
        type: File?
        inputBinding:
            position: 7
            prefix: --file-list
        doc: Supply a file containing a list of file paths to be searched, one per row.

outputs:
    outdir: 
        type: Directory?
        outputBinding:
            glob: $(inputs.filename + "_data")
    report:
        type: File
        outputBinding:
            glob: $(inputs.filename + ".html")
