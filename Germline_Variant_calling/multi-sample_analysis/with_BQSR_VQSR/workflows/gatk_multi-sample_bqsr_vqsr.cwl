class: Workflow

cwlVersion: v1.0

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement

inputs:
    ### Directory with fastq(.gz) input files ###
    raw_files_directory:
        type: Directory
    ### FASTQ file split ###
    input_file_split:
        type: string?
        default: "_R"
    input_file_split_fwd_single:
        type: string?
        default: "R1"
    input_file_split_rev:
        type: string?
        default: "R2"
    ### QC and trimming options ###
    input_qc_check: 
        type: boolean?
        default: true
    input_trimming_check:
        type: boolean?
        default: true
    ### Trim galore inputs ###
    tg_quality: int
    tg_length: int
    tg_compression: boolean
    tg_do_not_compress: boolean
    tg_strigency: int
    tg_trim_suffix: string
    ### Reference genome inputs ###
    reference_genome:
        type: File
        secondaryFiles:
            - .amb
            - .ann
            - .bwt
            - .pac
            - .sa
            - .fai
            - ^.dict
    ### BWA MEM inputs ###
    bwa_mem_sec_shorter_split_hits:
        type: boolean
        default: true
    bwa_mem_num_threads:
        type: int
        default: 16
    ### SAMtools inputs ###
    samtools_view_uncompressed:
        type: boolean
        default: false
    samtools_view_collapsecigar:
        type: boolean
        default: false
    samtools_view_readswithbits:
        type: int?
    samtools_view_readswithoutbits:
        type: int
        default: 0 # 2308
    samtools_view_fastcompression:
        type: boolean
        default: false
    samtools_view_samheader:
        type: boolean
        default: false
    samtools_view_count:
        type: boolean
        default: false
    samtools_view_readsingroup:
        type: string?
    samtools_view_readtagtostrip:
        type: string[]?
    samtools_view_readsquality:
        type: int?
    samtools_view_cigar:
        type: int?
    samtools_view_iscram:
        type: boolean
        default: false
    samtools_view_threads:
        type: int?
    samtools_view_randomseed:
        type: float?
    samtools_view_region:
        type: string?
    samtools_view_readsinlibrary:
        type: string?
    samtools_view_target_bed_file: 
        type: File?
    samtools_fixmate_threads:
        type: int
        default: 16
    samtools_fixmate_output_format:
        type: string
        default: bam
    samtools_sort_compression_level:
        type: int?
    samtools_sort_threads:
        type: int?
        default: 16
    samtools_sort_memory:
        type: string?
    samtools_flagstat_threads:
        type: int?
        default: 16
    ### Picard AddOrReplaceReadGroups inputs ###
    # picard_addorreplacereadgroups_rglb:
    #     type: string?
    picard_addorreplacereadgroups_rgpl:
        type: string?
    ### GATK - SplitIntervals inputs ###
    gatk_splitintervals_include_intervalList:
        type: File?
    gatk_splitintervals_exclude_intervalList: 
        type: File?
    gatk_splitintervals_scatter_count:
        type: int
    ### GATK - BQSR inputs ###
    # files associated to known variations
    sub_bqsr_known_sites_1:
        type: File
        secondaryFiles: 
            - .tbi
    sub_bqsr_known_sites_2:
        type: File?
        secondaryFiles: 
            - .tbi
    sub_bqsr_known_sites_3:
        type: File?
        secondaryFiles: 
            - .tbi
    sub_bqsr_interval_padding:
        type: int?
    sub_bqsr_hc_native_pairHMM_threads:
        type: int?
    sub_bqsr_hc_java_options: 
        type: string?
    ### GATK - VariantRecalibrator inputs ###
    VariantRecalibrator_use_annotation: 
        type: string[]
    VariantRecalibrator_trust_all_polymorphic:
        type: boolean?
    # INDELs
    VariantRecalibrator_truth_sensitivity_trance_indels:
        type: float[]?
    vqsr_arguments_indels_1:
        type: string
    vqsr_known_sites_indels_1:
        type: File
        secondaryFiles: 
            - .tbi
    vqsr_arguments_indels_2:
        type: string?
    vqsr_known_sites_indels_2:
        type: File?
        secondaryFiles: 
            - .tbi
    vqsr_arguments_indels_3:
        type: string?
    vqsr_known_sites_indels_3:
        type: File?
        secondaryFiles: 
            - .tbi
    # SNPs
    VariantRecalibrator_truth_sensitivity_trance_snps:
        type: float[]?
    vqsr_arguments_snps_1:
        type: string
    vqsr_known_sites_snps_1:
        type: File
        secondaryFiles: 
            - .tbi
    vqsr_arguments_snps_2:
        type: string?
    vqsr_known_sites_snps_2:
        type: File?
        secondaryFiles: 
            - .tbi
    vqsr_arguments_snps_3:
        type: string?
    vqsr_known_sites_snps_3:
        type: File?
        secondaryFiles: 
            - .tbi
    vqsr_arguments_snps_4:
        type: string?
    vqsr_known_sites_snps_4:
        type: File?
        secondaryFiles: 
            - .tbi
    ### GATK - ApplyVQSR inputs ###
    ApplyVQSR_ts_filter_level:
        type: float?
    ### bcftools view inputs ###
    bcftools_view_include_VQSR_filters:
        type: string
    bcftools_view_threads:
        type: int?
    ### bcftools norm inputs ###
    bcftools_norm_threads:
        type: int?
    bcftools_norm_multiallelics: 
        type: string
        default: "-both"
    ### ANNOVAR - table_annovar.pl inputs ###
    table_annovar_database_location:
        type: Directory
    table_annovar_build_over:
        type: string
    table_annovar_remove:
        type: boolean?
    table_annovar_protocol:
        type: string
    table_annovar_operation:
        type: string
    table_annovar_na_string:
        type: string?
    table_annovar_vcfinput:
        type: boolean
    table_annovar_otherinfo:
        type: boolean?
    table_annovar_convert_arg:
        type: string?

outputs:
    ### Trim_galore outputs ###
    o_trim_galore_single_fq:
        type: File[]
        outputSource: trim_galore_single/trim_galore
    o_trim_galore_single_reports:
        type: File[]
        outputSource: trim_galore_single/trim_galore_report
    o_trim_galore_paired_fq:
        type: File[]
        outputSource: trim_galore_paired/trim_galore
    o_trim_galore_paired_reports:
        type: File[]
        outputSource: trim_galore_paired/trim_galore_report
    ### FASTQC outputs ###
    # HTML
    o_fastqc_raw_html:
        type: File[]?
        outputSource: rename_fastqc_raw_html/renamed_file
    o_fastqc_single_html:
        type: File[]?
        outputSource: rename_fastqc_single_html/renamed_file
    o_fastqc_paired_html:
        type: File[]?
        outputSource: rename_fastqc_paired_html/renamed_file
    # ZIP
    o_fastqc_raw_zip:
        type: Directory?
        outputSource: cp_fastqc_raw_zip/output_dir
    o_fastqc_single_zip:
        type: Directory?
        outputSource: cp_fastqc_single_zip/output_dir
    o_fastqc_paired_zip:
        type: Directory?
        outputSource: cp_fastqc_paired_zip/output_dir
    ### BWA MEM outputs ###
    o_bwa_mem_single:
        type: File[]
        outputSource: bwa_mem_single/output
    o_bwa_mem_paired:
        type: File[]
        outputSource: bwa_mem_paired/output
    o_gather_bwa_sam_files:
        type: File[]
        outputSource: gather_bwa_sam_files/total_sam_files
    ### SAMtools outputs ###
    o_samtools_view_conversion:
        type: File[]
        outputSource: samtools_view_conversion/output
    samtools_sort_by_name:
        type: File[]
        outputSource: samtools_sort_by_name/sorted
    o_samtools_fixmate:
        type: File[]
        outputSource: samtools_fixmate/output
    o_samtools_sort:
        type: File[]
        outputSource: samtools_sort/sorted
    ### Picard tool outputs ###
    o_picard_addorreplacereadgroups:
        type: File[]
        outputSource: picard_addorreplacereadgroups/output
    o_picard_markduplicates:
        type: File[]
        outputSource: picard_markduplicates/output
    o_picard_markduplicates_metrics:
        type: File[]
        outputSource: picard_markduplicates/output_report
    ### SAMtools outputs ###
    o_samtools_flagstat:
        type: File[]
        outputSource: samtools_flagstat/output
    o_samtools_view_count_total:
        type: File[]
        outputSource: samtools_view_count_total/output
    o_samtools_index:
        type: File[]
        outputSource: samtools_index/alignments_with_index
    ### GATK SplitIntervals outputs ###
    o_gatk_splitintervals:
        type: File[]
        outputSource: gatk_splitintervals/intervalfiles
    ### GATK - BQSR outputs ###
    o_gatk_bqsr_subworkflowbqsr_tables:
        type: {
            "type": "array", 
            "items":{"type": "array", "items": "File"}
            }
        outputSource: gatk_bqsr_subworkflow/o_gatk_baserecalibrator_table
    o_gatk_bqsr_subworkflowbqsr_bqsr_bam:
        type: {
            "type": "array", 
            "items":{"type": "array", "items": "File"}
            }
        outputSource: gatk_bqsr_subworkflow/o_gatk_applybqsr
    o_gatk_bqsr_subworkflowbqsr_hc:
        type: {
            "type": "array", 
            "items":{"type": "array", "items": "File"}
            }
        outputSource: gatk_bqsr_subworkflow/o_gatk_HaplotypeCaller
    o_gatk_bqsr_subworkflowbqsr_mergevcfs:
        type: File[]
        outputSource: gatk_bqsr_subworkflow/o_gatk_MergeVCFs
    ### GATK - CombineGVCFs outputs ###
    o_gatk_CombineGVCFs:
        type: File
        outputSource: gatk_CombineGVCFs/output
    ### GATK - GenotypeGVCFs outputs ###
    o_gatk_GenotypeGVCFs:
        type: File
        outputSource: gatk_GenotypeGVCFs/output
    ### GATK - MakeSitesOnlyVcf ###
    o_gatk_MakeSitesOnlyVcf:
        type: File
        outputSource: gatk_MakeSitesOnlyVcf/output
    ### GATK - VariantRecalibrator ###
    o_gatk_VariantRecalibrator_indel_recal:
        type: File[]
        outputSource: gatk_VariantRecalibrator_indel/recal_table
    o_gatk_VariantRecalibrator_indel_tranches:
        type: File[]
        outputSource: gatk_VariantRecalibrator_indel/tranches_file
    o_gatk_VariantRecalibrator_snp_recal:
        type: File[]
        outputSource: gatk_VariantRecalibrator_snp/recal_table
    o_gatk_VariantRecalibrator_snp_tranches:
        type: File[]
        outputSource: gatk_VariantRecalibrator_snp/tranches_file
    ### GATK - ApplyVQSR ###
    o_gatk_ApplyVQSR_indel:
        type: File[]
        outputSource: gatk_ApplyVQSR_indel/output
    o_gatk_ApplyVQSR_snp:
        type: File[]
        outputSource: gatk_ApplyVQSR_snp/output
    
    ### MergeVCFs ###
    o_gatk_VQSR_MergeVCFs:
        type: File
        outputSource: gatk_VQSR_MergeVCFs/output
    # ### VQSR filtered VCF files - PASS only ###
    o_bcftools_view_filter_vqsr:
        type: File
        outputSource: bcftools_view_filter_vqsr/output
    ### bcftools norm for INDELs ###
    o_bcftools_norm_vqsr:
        type: File
        outputSource: bcftools_norm_vqsr/output   
    ### ANNOVAR outputs ###
    o_table_annovar_filtered_multianno_vcf:
        type: File
        outputSource: table_annovar_filtered/multianno_vcf
    o_table_annovar_filtered_multianno_txt:
        type: File
        outputSource: table_annovar_filtered/multianno_txt
    o_table_annovar_filtered_avinput:
        type: File
        outputSource: table_annovar_filtered/avinput

steps:
    ### Separate files - Generate file names (ExpressionTool) ###
    get_raw_files:
        run: ../wrappers/get-raw-files.cwl
        in:
            DIRECTORY: raw_files_directory
        out: [raw_files]
    ### Separate files - Generate file names (ExpressionTool) ###
    split_single_paired:
        run: ../wrappers/split-single-paired_v2.cwl
        in:
            qc_check: input_qc_check
            trimming_check: input_trimming_check
            input_files: get_raw_files/raw_files
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            file_split_rev: input_file_split_rev
        out: 
            - single_files
            - paired_files
            - trim_galore_for_single
            - trim_galore_for_paired
            - fastqc_for_raw
            - fastqc_for_single
            - fastqc_for_paired
            - cp_command_raw
            - cp_command_single
            - cp_command_paired
    ### Trim_galore ###
    trim_galore_single:
        run: ../wrappers/trim-galore.cwl
        in:
            command: split_single_paired/trim_galore_for_single
            fq_files: split_single_paired/single_files
            length: tg_length
            quality: tg_quality
            compression: tg_compression
            do_not_compress: tg_do_not_compress
            trim_suffix: tg_trim_suffix
            strigency: tg_strigency
            paired: 
                valueFrom: $( false )
        out:
            - trim_galore
            - trim_galore_report
    trim_galore_paired:
        run: ../wrappers/trim-galore.cwl
        in:
            command: split_single_paired/trim_galore_for_paired
            fq_files: split_single_paired/paired_files
            length: tg_length
            quality: tg_quality
            compression: tg_compression
            do_not_compress: tg_do_not_compress
            trim_suffix: tg_trim_suffix
            strigency: tg_strigency
            paired:
                valueFrom: $( true )
        out: 
            - trim_galore
            - trim_galore_report
    ### FASTQC ###
    fastqc_raw:
        run: ../wrappers/fastqc.cwl
        in:
            command: split_single_paired/fastqc_for_raw
            input_files: get_raw_files/raw_files
        out: 
            - html_file
            - zipped_file 
    fastqc_single_trimmed:
        run: ../wrappers/fastqc.cwl
        in:
            command: split_single_paired/fastqc_for_single
            input_files: trim_galore_single/trim_galore
        out:
            - html_file
            - zipped_file
    fastqc_paired_trimmed:
        run: ../wrappers/fastqc.cwl
        in:
            command: split_single_paired/fastqc_for_paired
            input_files: trim_galore_paired/trim_galore
        out:
            - html_file
            - zipped_file
    ### Copy FASTQC output files to directory - ZIP ###
    cp_fastqc_raw_zip:
        run: ../wrappers/cp.cwl
        in:
            command: split_single_paired/cp_command_raw
            input_files: fastqc_raw/zipped_file
            outputdir:
                valueFrom: $( "fastqc_raw_zip" )
        out:
            - output_dir
    cp_fastqc_single_zip:
        run: ../wrappers/cp.cwl
        in:
            command: split_single_paired/cp_command_single
            input_files: fastqc_single_trimmed/zipped_file
            outputdir:
                valueFrom: $( "fastqc_single_trimmed_zip" )
        out:
            - output_dir
    cp_fastqc_paired_zip:
        run: ../wrappers/cp.cwl
        in:
            command: split_single_paired/cp_command_paired
            input_files: fastqc_paired_trimmed/zipped_file
            outputdir:
                valueFrom: $( "fastqc_paired_trimmed_zip" )
        out:
            - output_dir
    ### Rename FASTQC output files - HTML ###
    rename_fastqc_raw_html:
        run: ../wrappers/rename.cwl
        in:
            input_file: fastqc_raw/html_file
            run_type:
                valueFrom: $( "fastqc_raw_" )
        out:
            - renamed_file
        scatter:
            - input_file
    rename_fastqc_single_html:
        run: ../wrappers/rename.cwl
        in:
            input_file: fastqc_single_trimmed/html_file
            run_type:
                valueFrom: $( "fastqc_single_trimmed_" )
        out:
            - renamed_file
        scatter:
            - input_file
    rename_fastqc_paired_html:
        run: ../wrappers/rename.cwl
        in:
            input_file: fastqc_paired_trimmed/html_file
            run_type:
                valueFrom: $( "fastqc_paired_trimmed_" )
        out:
            - renamed_file
        scatter:
            - input_file    
    ### Check for trimming option ###
    check_trimming:
        run: ../wrappers/check_trimming.cwl
        in: 
            trimming_check: input_trimming_check
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            file_split_rev: input_file_split_rev
            input_single: split_single_paired/single_files
            input_paired: split_single_paired/paired_files
            trimming_single: trim_galore_single/trim_galore
            trimming_paired: trim_galore_paired/trim_galore
        out:
        - single_files
        - paired_files
    ### Collect flowcell ID and lane information - Single end ###
    rg_extraction_single:
        run: ../wrappers/fastq_RG_extraction.cwl
        scatter:
            - input
        in: 
            input: check_trimming/single_files
        out: [rg_information]
    ### bwa mem - Single end ###
    bwa_mem_single:
        run: ../wrappers/bwa-mem.cwl
        scatter:
        - trimmed_fq_read1
        in:
            file_split: input_file_split
            sec_shorter_split_hits: bwa_mem_sec_shorter_split_hits
            num_threads: bwa_mem_num_threads
            ref_genome: reference_genome
            trimmed_fq_read1: check_trimming/single_files
        out: [output]
    ### Split forward and reverse read files ###
    split_paired_read1_read2:
        run: ../wrappers/split-paired-read1-read2.cwl
        in:
            paired_files: check_trimming/paired_files
            file_split_fwd_single: input_file_split_fwd_single
            file_split_rev: input_file_split_rev
        out: [reads_1, reads_2]
    ### Collect flowcell ID and lane information - Paired end ###
    rg_extraction_paired:
        run: ../wrappers/fastq_RG_extraction.cwl
        scatter:
            - input
        in: 
            input: split_paired_read1_read2/reads_1
        out: [rg_information]
    ### bwa mem - Paired end ###
    bwa_mem_paired:
        run: ../wrappers/bwa-mem.cwl
        scatterMethod: dotproduct
        scatter:
        - trimmed_fq_read1
        - trimmed_fq_read2
        in:
            file_split: input_file_split
            sec_shorter_split_hits: bwa_mem_sec_shorter_split_hits
            num_threads: bwa_mem_num_threads
            ref_genome: reference_genome
            trimmed_fq_read1: split_paired_read1_read2/reads_1
            trimmed_fq_read2: split_paired_read1_read2/reads_2
        out: [output]
    ### Collect and output all SAM files & generate filenames for next steps ###
    gather_bwa_sam_files:
        run: ../wrappers/gather-bwa-sam-files.cwl
        in: 
            single_files: bwa_mem_single/output
            paired_files: bwa_mem_paired/output
            single_files_rg_info: rg_extraction_single/rg_information
            paired_files_rg_info: rg_extraction_paired/rg_information
        out: [total_sam_files,
              names_rgids,
              names_rgpus,
              names_rgplbs]
    ### samtools view - SAM to BAM conversion ###
    # options for filtering out based on SAM FLAG fields:
    # read unmapped (0x4)
    # not primary alignment (0x100)
    # supplementary alignment (0x800)
    samtools_view_conversion:
        run: ../wrappers/samtools-view.cwl
        scatter:
        - input
        in:
            input: gather_bwa_sam_files/total_sam_files
            isbam:
                valueFrom: $( true )
            collapsecigar: samtools_view_collapsecigar
            readsingroup: samtools_view_readsingroup
            uncompressed: samtools_view_uncompressed
            readtagtostrip: samtools_view_readtagtostrip
            readsquality: samtools_view_readsquality
            ### 
            # keep or remove reads based on SAM FLAG values:
            readswithbits: samtools_view_readswithbits
            readswithoutbits: samtools_view_readswithoutbits
            ###
            cigar: samtools_view_cigar
            iscram: samtools_view_iscram
            threads: samtools_view_threads
            fastcompression: samtools_view_fastcompression
            samheader: samtools_view_samheader
            count: samtools_view_count
            randomseed: samtools_view_randomseed
            region: samtools_view_region
            readsinlibrary: samtools_view_readsinlibrary
            output_name: 
                valueFrom: $( inputs.input.basename.split(".sam")[0].concat(".bam") )
            # (optional) filter also with genomics coordinates:
            target_bed_file: samtools_view_target_bed_file
        out: [output]
    ### samtools sort by read names ###
    samtools_sort_by_name:
        run: ../wrappers/samtools-sort.cwl
        scatter:
        - input
        in:
            compression_level: samtools_sort_compression_level 
            threads: samtools_sort_threads 
            memory: samtools_sort_memory 
            input: samtools_view_conversion/output
            output_name: 
                valueFrom: $( inputs.input.basename.split(".bam")[0].concat(".name.sorted.bam") )
            sort_by_name:
             valueFrom: $( true )
        out: [sorted]
    ### samtools fixmate – fills in mate coordinates and insert size fields ###
    samtools_fixmate:
        run: ../wrappers/samtools-fixmate.cwl
        scatter:
        - input_file
        in:
            threads: samtools_fixmate_threads
            output_format: samtools_fixmate_output_format
            input_file: samtools_sort_by_name/sorted
            output_file_name: 
                valueFrom: $( inputs.input_file.basename.split(".name.sorted.bam")[0].concat(".fixed.bam") )
        out: [output]
    ### samtools sort by chromosomal coordinates ###
    samtools_sort:
        run: ../wrappers/samtools-sort.cwl
        scatter:
        - input
        in:
            compression_level: samtools_sort_compression_level 
            threads: samtools_sort_threads 
            memory: samtools_sort_memory 
            input: samtools_fixmate/output
            output_name: 
                valueFrom: $( inputs.input.basename.split(".fixed.bam")[0].concat(".fixed.sorted.bam") )
        out: [sorted]
    ### Picard AddOrReplaceReadGroups ###
    picard_addorreplacereadgroups:
        run: ../wrappers/picard-AddOrReplaceReadGroups.cwl
        scatterMethod: dotproduct
        scatter:
        - INPUT
        - rgid
        - rgpu
        - rglb
        in: 
            INPUT: samtools_sort/sorted 
            OUTPUT: 
                valueFrom: $( inputs.INPUT.basename.split(".fixed.sorted.bam")[0].concat(".fixed.sorted.uniq.rg.bam") )
            rgid: gather_bwa_sam_files/names_rgids
            rglb: gather_bwa_sam_files/names_rgplbs
            rgpl: picard_addorreplacereadgroups_rgpl 
            rgpu: gather_bwa_sam_files/names_rgpus
            rgsm: 
                valueFrom: $( inputs.INPUT.basename.split(".fixed.sorted.bam")[0] )
        out: [output]
    ### Picard MarkDuplicates ###
    picard_markduplicates:
        run: ../wrappers/picard-MarkDuplicates.cwl
        scatter: 
        - INPUT
        in:
            INPUT: picard_addorreplacereadgroups/output
            OUTPUT: 
                valueFrom: $( inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.bam")[0].concat(".fixed.sorted.uniq.rg.md.bam") )
            METRICS_FILE: 
                valueFrom: $( inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.bam")[0].concat(".md_metrics.txt") )
        out: [output, output_report]
    ### samtools flagstat – counts the number of alignments for each FLAG type ###
    samtools_flagstat:
        run: ../wrappers/samtools-flagstat.cwl
        scatter:
        - input
        in:
            threads: samtools_flagstat_threads
            input: picard_markduplicates/output
            output_name: 
                valueFrom: $( inputs.input.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".align.stats.txt") )
        out: [output]
    ### samtools view - print only alignment counts ###
    samtools_view_count_total:
        run: ../wrappers/samtools-view.cwl
        scatter:
        - input
        in:
            isbam: 
                valueFrom: $( false )
            collapsecigar: samtools_view_collapsecigar
            readsingroup: samtools_view_readsingroup
            uncompressed: samtools_view_uncompressed
            readtagtostrip: samtools_view_readtagtostrip
            input: picard_markduplicates/output
            readsquality: samtools_view_readsquality
            # readswithbits: samtools_view_readswithbits
            # readswithoutbits: samtools_view_readswithoutbits
            cigar: samtools_view_cigar
            iscram: samtools_view_iscram
            threads: samtools_view_threads
            fastcompression: samtools_view_fastcompression
            samheader: samtools_view_samheader
            count: 
                default: true
            randomseed: samtools_view_randomseed
            region: samtools_view_region
            readsinlibrary: samtools_view_readsinlibrary
            output_name: 
                valueFrom: $( inputs.input.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".count.total.txt") )
            # target_bed_file: samtools_view_target_bed_file
        out: [output]
    ### GATK SplitIntervals ###
    gatk_splitintervals:
        run: ../wrappers/gatk-SplitIntervals.cwl
        in:
            reference: reference_genome
            include_intervalList: gatk_splitintervals_include_intervalList
            exclude_intervalList: gatk_splitintervals_exclude_intervalList
            scatter_count: gatk_splitintervals_scatter_count
        out: [intervalfiles]
    ### samtools index ###
    samtools_index:
        run: ../wrappers/samtools-index.cwl
        scatter:
        - alignments
        in:
            alignments: picard_markduplicates/output
        out: [alignments_with_index]
    ### GATK BaseRecalibrator - subworkflow ###
    gatk_bqsr_subworkflow:
        run: ../workflows/gatk-bqsr-subworkflow_multiple.cwl
        # scatterMethod: dotproduct
        scatter: 
        - bqsr_INPUT
        # - bqsr_OUTPUT
        in: 
            bqsr_hc_java_options: sub_bqsr_hc_java_options 
            bqsr_INPUT: samtools_index/alignments_with_index 
            bqsr_OUTPUT: 
                valueFrom: $( inputs.bqsr_INPUT.basename.split(".fixed.sorted.uniq.rg.md.bam")[0] )
            bqsr_reference: reference_genome
            bqsr_known_sites_1: sub_bqsr_known_sites_1
            bqsr_known_sites_2: sub_bqsr_known_sites_2
            bqsr_known_sites_3: sub_bqsr_known_sites_3
            bqsr_intervals: gatk_splitintervals/intervalfiles
            bqsr_interval_padding: sub_bqsr_interval_padding
            bqsr_hc_native_pairHMM_threads: sub_bqsr_hc_native_pairHMM_threads
        out: 
        - o_gatk_baserecalibrator_table
        - o_gatk_applybqsr
        - o_gatk_HaplotypeCaller
        - o_gatk_MergeVCFs
    ### GATK - CombineGVCFs ###
    gatk_CombineGVCFs:
        run: ../wrappers/gatk-CombineGVCFs.cwl
        in: 
            reference: reference_genome
            gvcf_files: gatk_bqsr_subworkflow/o_gatk_MergeVCFs
            output_name:
                valueFrom: $( "cohort.g.vcf.gz" )
        out: [output]
    ### GATK - GenotypeGVCFs ###
    gatk_GenotypeGVCFs:
        run: ../wrappers/gatk-GenotypeGVCFs.cwl 
        in:
            reference: reference_genome
            gvcf_input: gatk_CombineGVCFs/output
            vcf_output: 
                valueFrom: $( "cohort.genotyped.vcf.gz" )
        out: [output]
    # GATK - VQSR
    ### GATK - MakeSitesOnlyVcf ###
    gatk_MakeSitesOnlyVcf:
        run: ../wrappers/gatk-MakeSitesOnlyVcf.cwl 
        in: 
            INPUT: gatk_GenotypeGVCFs/output
            OUTPUT: 
                valueFrom: $( "cohort.genotyped.sitesonly.vcf.gz" )
        out: [output]
    ### GATK - VariantRecalibrator - INDELs ###
    gatk_VariantRecalibrator_indel:
        run: ../wrappers/gatk-VariantRecalibrator.cwl
        scatter: 
        - intervals
        in: 
            truth_sensitivity_trance: VariantRecalibrator_truth_sensitivity_trance_indels
            use_annotation: VariantRecalibrator_use_annotation
            reference: reference_genome
            variant: gatk_MakeSitesOnlyVcf/output
            mode: 
                valueFrom: $( "INDEL" )
            arguments_1: vqsr_arguments_indels_1
            resource_1: vqsr_known_sites_indels_1
            arguments_2: vqsr_arguments_indels_2
            resource_2: vqsr_known_sites_indels_2
            arguments_3: vqsr_arguments_indels_3
            resource_3: vqsr_known_sites_indels_3
            intervals: gatk_splitintervals/intervalfiles
            trust_all_polymorphic: VariantRecalibrator_trust_all_polymorphic
            output_name: 
                valueFrom: $( "cohort.genotyped.indel" )
        out:
        - recal_table
        - tranches_file
    ### GATK - VariantRecalibrator - SNPs ###
    gatk_VariantRecalibrator_snp:
        run: ../wrappers/gatk-VariantRecalibrator.cwl
        scatter: 
        - intervals
        in: 
            truth_sensitivity_trance: VariantRecalibrator_truth_sensitivity_trance_snps
            use_annotation: VariantRecalibrator_use_annotation
            reference: reference_genome
            variant: gatk_MakeSitesOnlyVcf/output
            mode: 
                valueFrom: $( "SNP" )
            arguments_1: vqsr_arguments_snps_1
            resource_1: vqsr_known_sites_snps_1
            arguments_2: vqsr_arguments_snps_2
            resource_2: vqsr_known_sites_snps_2
            arguments_3: vqsr_arguments_snps_3
            resource_3: vqsr_known_sites_snps_3
            arguments_4: vqsr_arguments_snps_4
            resource_4: vqsr_known_sites_snps_4
            intervals: gatk_splitintervals/intervalfiles
            trust_all_polymorphic: VariantRecalibrator_trust_all_polymorphic
            output_name: 
                valueFrom: $( "cohort.genotyped.snp" )
        out:
        - recal_table
        - tranches_file
    ### GATK - ApplyVQSR - INDELs ###
    gatk_ApplyVQSR_indel: 
        run: ../wrappers/gatk-ApplyVQSR.cwl
        scatterMethod: dotproduct
        scatter:
        - recal_file
        - tranches_file
        - intervals
        in:
            reference: reference_genome
            variant: gatk_GenotypeGVCFs/output
            ts_filter_level: ApplyVQSR_ts_filter_level
            recal_file: gatk_VariantRecalibrator_indel/recal_table
            tranches_file: gatk_VariantRecalibrator_indel/tranches_file
            intervals: gatk_splitintervals/intervalfiles
            mode:
                valueFrom: $( "INDEL" )
            output_name: 
                valueFrom: $( "indel.recalibrated" )
        out: [output]
    ### GATK - ApplyVQSR - SNPs ###
    gatk_ApplyVQSR_snp: 
        run: ../wrappers/gatk-ApplyVQSR.cwl
        scatterMethod: dotproduct
        scatter:
        - variant
        - recal_file
        - tranches_file
        - intervals
        in:
            reference: reference_genome
            variant: gatk_ApplyVQSR_indel/output
            ts_filter_level: ApplyVQSR_ts_filter_level
            recal_file: gatk_VariantRecalibrator_snp/recal_table
            tranches_file: gatk_VariantRecalibrator_snp/tranches_file
            intervals: gatk_splitintervals/intervalfiles
            mode:
                valueFrom: $( "SNP" )
            output_name: 
                valueFrom: $( "snp.recalibrated" )
        out: [output]
    gatk_VQSR_MergeVCFs:
        run: ../wrappers/gatk-MergeVCFs_vqsr.cwl
        in:
            INPUT: gatk_ApplyVQSR_snp/output
            OUTPUT: 
                valueFrom: $( "cohort.snp.indel.vqsr.recalibrated" )
        out: [output]
    ### bcftools view - filter variants ###
    bcftools_view_filter_vqsr:
        run: ../wrappers/bcftools-view.cwl
        in: 
            input: gatk_VQSR_MergeVCFs/output
            threads: bcftools_view_threads
            include: bcftools_view_include_VQSR_filters
            output_type: 
                valueFrom: $( "z" )
            output_name:
                valueFrom: $("cohort.snp.indel.vqsr.recalibrated.filtered.vcf.gz")
        out: [output]
    ### bcftools norm ###
    bcftools_norm_vqsr:
        run: ../wrappers/bcftools-norm.cwl
        in: 
            input: gatk_VQSR_MergeVCFs/output
            threads: bcftools_norm_threads
            reference: reference_genome
            multiallelics: bcftools_norm_multiallelics
            output_type: 
                valueFrom: $( "v" )
        out: [output]
    ### ANNOVAR - variant annotation ###
    table_annovar_filtered:
        run: ../wrappers/table-annovar.cwl
        in:
            query_file: bcftools_norm_vqsr/output
            database_location: table_annovar_database_location
            build_over: table_annovar_build_over
            remove: table_annovar_remove
            protocol: table_annovar_protocol
            operation: table_annovar_operation
            na_string: table_annovar_na_string
            vcfinput: table_annovar_vcfinput
            otherinfo: table_annovar_otherinfo
            convert_arg: table_annovar_convert_arg
        out:
        - multianno_vcf
        - multianno_txt
        - avinput
