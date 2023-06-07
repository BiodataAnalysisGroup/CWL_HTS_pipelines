cwlVersion: v1.0
class: CommandLineTool

requirements:
    - class: DockerRequirement
      dockerPull: ubuntu:latest
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement

arguments:
    - zcat
    - $(inputs.input)
    - shellQuote: false
      valueFrom: '|'
    - grep
    - "@"
    - shellQuote: false
      valueFrom: '|'
    - head
    - "-1"
    - shellQuote: false
      valueFrom: '|'
    - cut
    - shellQuote: false
      valueFrom: '-d ":" -f 3-4'
    - shellQuote: false
      valueFrom: '|'
    - sed
    - 's/:/./'
    - shellQuote: false
      valueFrom: '|'
    - awk
    - '{print "$(inputs.input.basename.split(inputs.file_split)[0])."$1}'
    - shellQuote: false
      valueFrom: '>'
    # - shellQuote: false
    #   valueFrom: '|'
    # - tee
    # - "-a"
    - $(inputs.input.basename.split(inputs.file_split)[0] + "_rg_info.txt")
    # - $(inputs.output_file)

inputs: 
    input: 
        type: File
    file_split:
        type: string
        default: "_R"
    # output_file:
    #     type: string?
    #     default: "rg_info.txt"

outputs:
    rg_information: 
        type: string #File
        outputBinding:
          # glob: $(inputs.output_file)
          glob: $(inputs.input.basename.split(inputs.file_split)[0] + "_rg_info.txt")
          loadContents: true
          outputEval: $(self[0].contents.split("\n")[0])
