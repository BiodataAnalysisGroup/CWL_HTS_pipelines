cwlVersion: v1.0
class: ExpressionTool

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: single_files
  type: File[]
- id: paired_files
  type: File[]

outputs:
- id: total_sam_files
  type: File[]
- id: names_basenames
  type: string[]
- id: names_bam_raw
  type: string[]
- id: names_bam_fixed
  type: string[]
- id: names_bam_sorted
  type: string[]
- id: names_bam_uniq
  type: string[]
- id: names_bam_uniq_rg
  type: string[]
- id: names_txt_align_stats
  type: string[]
- id: names_txt_coverage_mean
  type: string[]
- id: names_txt_count_ontarget
  type: string[]
- id: names_txt_count_total
  type: string[]
- id: names_rgids
  type: int[]
- id: names_rgpus
  type: string[]

expression: |
  ${
    var out = [];

    var bam_raw = [];
    var bam_fixed = [];
    var bam_sorted = [];
    var bam_uniq = [];
    var bam_uniq_rg = [];

    var txt_align_stats = [];
    var txt_coverage_mean = []
    var txt_count_ontarget = [];
    var txt_count_total = [];

    var picard_rgids = [];
    var picard_rgpus = [];

    var basenames = [];

    var i;
    var sample;
    var num;

    for(i = 0; i < inputs.single_files.length; i++){
      out.push(inputs.single_files[i]);
    }

    for(i = 0; i < inputs.paired_files.length; i++){
      out.push(inputs.paired_files[i]);
    }

    for(i = 0; i < out.length; i++){
      sample = out[i].basename.split(".sam")[0];

      bam_raw.push(sample.concat(".bam"));
      bam_fixed.push(sample.concat(".fixed.bam"));
      bam_sorted.push(sample.concat(".fixed.sorted.bam"));
      bam_uniq.push(sample.concat(".fixed.sorted.uniq.bam"));
      bam_uniq_rg.push(sample.concat(".fixed.sorted.uniq.rg.bam"));
      txt_align_stats.push(sample.concat(".align.stats.txt"));
      txt_coverage_mean.push(sample.concat(".coverage.mean.txt"));
      txt_count_ontarget.push(sample.concat(".count.ontarget.txt"));
      txt_count_total.push(sample.concat(".count.total.txt"));
      basenames.push(sample);


      num = i + 1;
      picard_rgids.push(num);

      sample = num.toString();
      picard_rgpus.push("unit".concat(sample));
    }

    return {
      "total_sam_files": out,
      "names_basenames": basenames,
      "names_bam_raw": bam_raw,
      "names_bam_fixed": bam_fixed,
      "names_bam_sorted": bam_sorted,
      "names_bam_uniq": bam_uniq,
      "names_bam_uniq_rg": bam_uniq_rg,
      "names_txt_align_stats": txt_align_stats,
      "names_txt_coverage_mean": txt_coverage_mean,
      "names_txt_count_ontarget": txt_count_ontarget,
      "names_txt_count_total": txt_count_total,
      "names_rgids": picard_rgids,
      "names_rgpus": picard_rgpus
    };
  }