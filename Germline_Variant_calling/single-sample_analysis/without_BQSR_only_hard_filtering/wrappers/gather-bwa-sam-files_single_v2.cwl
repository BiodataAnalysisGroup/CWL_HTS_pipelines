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
- id: names_bam_uniq_rg_md
  type: string[]
- id: names_bam_uniq_rg_md_metrics
  type: string[]
# - id: names_bqsr_tables
#   type: string[]
- id: names_bam_uniq_rg_md_bqsr
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
  type: string[]
- id: names_rgpus
  type: string[]
- id: names_rgplbs
  type: string[]
- id: names_raw_vcf
  type: string[]
- id: names_snp_vcf
  type: string[]
- id: names_indel_vcf
  type: string[]
- id: names_filtered_snp_vcf
  type: string[]
- id: names_filtered_indel_vcf
  type: string[]
- id: names_bcftools_concat_vcf
  type: string[]
- id: names_bcftools_hard_filtered_vcf
  type: string[]
- id: names_cnn_vcf
  type: string[]
- id: names_cnn_filtered_vcf
  type: string[]
- id: names_bcftools_cnn_filtered_vcf
  type: string[]

expression: |
  ${
    var out = [];
    var out_rg = [];

    var bam_raw = [];
    var bam_fixed = [];
    var bam_sorted = [];
    var bam_uniq = [];
    var bam_uniq_rg = [];
    var bam_uniq_rg_md = [];
    var bam_uniq_rg_md_metrics = [];

    var bqsr_tables = [];
    var bam_uniq_rg_md_bqsr = [];

    var txt_align_stats = [];
    var txt_coverage_mean = []
    var txt_count_ontarget = [];
    var txt_count_total = [];

    var picard_rgids = [];
    var picard_rgpus = [];
    var picard_rglbs = [];

    var raw_vcf = [];
    var snp_vcf = [];
    var indel_vcf = [];
    var filtered_snp_vcf = [];
    var filtered_indel_vcf = [];
    var bcftools_concat_vcf = [];
    var bcftools_hard_filtered_vcf = [];

    var cnn_vcf = [];
    var cnn_filtered_vcf = [];
    var bcftools_cnn_filtered_vcf = [];

    var basenames = [];

    var i;
    var sample;
    var num;

    // SAM/BAM file information
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
      bam_uniq_rg.push(sample.concat(".fixed.sorted.uniq.rg.bam"));
      bam_uniq_rg_md.push(sample.concat(".fixed.sorted.uniq.rg.md.bam"));

      bam_uniq.push(sample.concat(".fixed.sorted.uniq.rg.md.filtered.bam")); // samtools view -> filter based on SAM FLAGs

      // bqsr_tables.push(sample.concat(".bqsr.table"));
      bam_uniq_rg_md_bqsr.push(sample.concat(".bqsr.bam"));

      bam_uniq_rg_md_metrics.push(sample.concat(".md_metrics.txt"));
      txt_align_stats.push(sample.concat(".align.stats.txt"));
      txt_coverage_mean.push(sample.concat(".coverage.mean.txt"));
      txt_count_ontarget.push(sample.concat(".count.ontarget.txt"));
      txt_count_total.push(sample.concat(".count.total.txt"));

      raw_vcf.push(sample.concat(".raw.vcf"));
      snp_vcf.push(sample.concat(".snp.vcf"));
      indel_vcf.push(sample.concat(".indel.vcf"));
      filtered_snp_vcf.push(sample.concat(".filtered.snp.vcf"));
      filtered_indel_vcf.push(sample.concat(".filtered.indel.vcf"));
      bcftools_concat_vcf.push(sample.concat(".concat.vcf"));
      bcftools_hard_filtered_vcf.push(sample.concat(".bcftools.hard.filtered.vcf"));
      
      cnn_vcf.push(sample.concat(".cnn.vcf"));
      cnn_filtered_vcf.push(sample.concat(".cnn_filtered.vcf"));
      bcftools_cnn_filtered_vcf.push(sample.concat(".bcftools_cnn_filtered.vcf"));
      
      basenames.push(sample);
    }

    // RG information
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

    return {
      "total_sam_files": out,
      "names_basenames": basenames,
      "names_bam_raw": bam_raw,
      "names_bam_fixed": bam_fixed,
      "names_bam_sorted": bam_sorted,
      "names_bam_uniq": bam_uniq,
      "names_bam_uniq_rg": bam_uniq_rg,
      "names_bam_uniq_rg_md": bam_uniq_rg_md,
      "names_bam_uniq_rg_md_metrics": bam_uniq_rg_md_metrics,
      // "names_bqsr_tables": bqsr_tables,
      "names_bam_uniq_rg_md_bqsr": bam_uniq_rg_md_bqsr,
      "names_txt_align_stats": txt_align_stats,
      "names_txt_coverage_mean": txt_coverage_mean,
      "names_txt_count_ontarget": txt_count_ontarget,
      "names_txt_count_total": txt_count_total,
      "names_rgids": picard_rgids,
      "names_rgpus": picard_rgpus,
      "names_rgplbs": picard_rglbs,
      "names_raw_vcf": raw_vcf,
      "names_snp_vcf": snp_vcf,
      "names_indel_vcf": indel_vcf,
      "names_filtered_snp_vcf": filtered_snp_vcf,
      "names_filtered_indel_vcf": filtered_indel_vcf,
      "names_bcftools_concat_vcf": bcftools_concat_vcf,
      "names_bcftools_hard_filtered_vcf": bcftools_hard_filtered_vcf,
      "names_cnn_vcf": cnn_vcf,
      "names_cnn_filtered_vcf": cnn_filtered_vcf,
      "names_bcftools_cnn_filtered_vcf": bcftools_cnn_filtered_vcf
    };
  }
