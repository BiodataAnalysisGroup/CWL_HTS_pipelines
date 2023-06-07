cwlVersion: v1.0
class: ExpressionTool

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: single_files
  type: File[]
- id: paired_files
  type: File[]
- id: single_files_rg_info
  type: string[]
- id: paired_files_rg_info
  type: string[]

outputs:
- id: total_sam_files
  type: File[]
- id: names_rgids
  type: string[]
- id: names_rgpus
  type: string[]
- id: names_rgplbs
  type: string[]

expression: |
  ${
    // set variables
    var out = [];
    // RG information variables
    var out_rg = [];
    var picard_rgids = [];
    var picard_rgpus = [];
    var picard_rglbs = [];
    // iteration variable
    var i;

    // SAM/BAM file information
    for(i = 0; i < inputs.single_files.length; i++){
      out.push(inputs.single_files[i]);
    }
    for(i = 0; i < inputs.paired_files.length; i++){
      out.push(inputs.paired_files[i]);
    }

    // Extract RG information for picard AddOrReplaceReadGroups
    for(i = 0; i < inputs.single_files_rg_info.length; i++){
      out_rg.push(inputs.single_files_rg_info[i]);
    }
    for(i = 0; i < inputs.paired_files_rg_info.length; i++){
      out_rg.push(inputs.paired_files_rg_info[i]);
    }
    for(i = 0; i < out_rg.length; i++){
      // RG IDs
      picard_rgids.push(out_rg[i]);
      // RG LBs
      picard_rglbs.push(out_rg[i].split(".")[0]);
      // RG PUs
      let flowcellID = out_rg[i].split(".")[1];
      let flowcellLane = out_rg[i].split(".")[2];
      picard_rgpus.push( flowcellID + "." + flowcellLane);
    }

    // return results
    return {
      "total_sam_files": out,
      "names_rgids": picard_rgids,
      "names_rgpus": picard_rgpus,
      "names_rgplbs": picard_rglbs
    };
  }
