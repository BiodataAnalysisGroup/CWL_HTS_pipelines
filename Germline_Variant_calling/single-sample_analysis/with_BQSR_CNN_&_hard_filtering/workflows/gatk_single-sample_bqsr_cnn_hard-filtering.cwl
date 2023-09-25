class: Workflow

cwlVersion: v1.0

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: SubworkflowFeatureRequirement

inputs:
    # Directory with fastq(.gz) input files
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
    samtools_view_readswithoutbits:
        type: int
        default: 0x904
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
    samtools_view_readswithbits:
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
    ### GATK - SplitIntervals ###
    gatk_splitintervals_include_intervalList:
        type: File?
    gatk_splitintervals_exclude_intervalList: 
        type: File?
    gatk_splitintervals_scatter_count:
        type: int
    ### GATK - BQSR inputs ###
    sub_bqsr_known_sites_1:
        type: File
        secondaryFiles: 
            - .tbi
    sub_bqsr_known_sites_2:
        type: File
        secondaryFiles: 
            - .tbi
    sub_bqsr_known_sites_3:
        type: File
        secondaryFiles: 
            - .tbi
    sub_bqsr_interval_padding:
        type: int?
    ### GATK - Apply-BQSR inputs ###
    ### GATK HaplotypeCaller inputs ###
    sub_hc_native_pairHMM_threads:
        type: int?  
    sub_hc_java_options:
        type: string?
    ### GATK SelectVariants inputs ###
    ### GATK VariantFiltration inputs ###
    VariantFiltration_window:
        type: int
        default: 0
    VariantFiltration_cluster:
        type: int
        default: 3
    # SNPs
    VariantFiltration_filter_name_snp:
        type:
            type: array
            items: string
    VariantFiltration_filter_snp:
        type:
            type: array
            items: string
    # INDELs
    VariantFiltration_filter_name_indel:
        type:
            type: array
            items: string
    VariantFiltration_filter_indel:
        type:
            type: array
            items: string
    ### GATK CNNScoreVariants inputs ###
    ### GATK FilterVariantTranches inputs ###
    FilterVariantTranches_resource_1:
        type: File
        secondaryFiles: 
        - .tbi
    FilterVariantTranches_resource_2:
        type: File?
        secondaryFiles: 
        - .tbi
    FilterVariantTranches_resource_3:
        type: File?
        secondaryFiles: 
        - .tbi
    ### bcftools view inputs ###
    bcftools_view_include_hard_filters:
        type: string
    bcftools_view_include_CNN_filters:
        type: string
    bcftools_view_threads:
        type: int
    ### bcftools norm inputs ###
    bcftools_norm_threads:
        type: int?
    bcftoomls_norm_multiallelics: 
        type: string
        default: "-both"
    ### ANNOVAR table_annovar.pl inputs ###
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
    ### SAMtools outputs ###
    o_gather_bwa_sam_files:
        type: File[]
        outputSource: gather_bwa_sam_files/total_sam_files
    o_samtools_view_conversion:
        type: File[]
        outputSource: samtools_view_conversion/output
    o_samtools_sort_by_name:
        type: File[]
        outputSource: samtools_sort_by_name/sorted
    o_samtools_fixmate:
        type: File[]
        outputSource: samtools_fixmate/output
    o_samtools_sort:
        type: File[]
        outputSource: samtools_sort/sorted
    o_picard_addorreplacereadgroups:
        type: File[]
        outputSource: picard_addorreplacereadgroups/output
    o_picard_markduplicates:
        type: File[]
        outputSource: picard_markduplicates/output
    o_picard_markduplicates_metrics:
        type: File[]
        outputSource: picard_markduplicates/output_report
    o_samtools_flagstat:
        type: File[]
        outputSource: samtools_flagstat/output
    o_samtools_view_count_total:
        type: File[]
        outputSource: samtools_view_count_total/output
    o_samtools_index:
        type: File[]
        outputSource: samtools_index/alignments_with_index
    o_gatk_bqsr_subworkflow:
        type: File[]
        outputSource: gatk_bqsr_subworkflow/o_gatk_gatherbqsrreports
    o_gatk_ApplyBQSR:
        type: File[]
        outputSource: gatk_applybqsr/output
    o_samtools_index_2:
        type: File[]
        outputSource: samtools_index_2/alignments_with_index
    # GATK SplitIntervals outputs
    o_gatk_splitintervals:
        type: File[]
        outputSource: gatk_splitintervals/intervalfiles
    ### VCF output files ###
    # Raw VCF files
    o_gatk_HaplotypeCaller:
        type: File[]
        outputSource: gatk_haplotypecaller_subworkflow/o_gatk_MergeVCFs
    # Hard-filtered SNP & INDEL VCF files
    o_tabix_snps:
        type: File[]
        outputSource: tabix_snps/output
    o_tabix_indels:
        type: File[]
        outputSource: tabix_indels/output
    # Concatenated, hard-filtered VCF files with FILTER tags
    o_bcftools_concat:
        type: File[]
        outputSource: bcftools_concat/output
    # Concatenated, hard-filtered VCF files - PASS only
    o_bcftools_view_hard_filter:
        type: File[]
        outputSource: bcftools_view_hard_filter/output
    o_bcftools_norm_hard_filter:
        type: File[]
        outputSource: bcftools_norm_hard_filter/output
    # CNNScoreVariants VCF files with scores
    o_gatk_CNNScoreVariants:
        type: File[]
        outputSource: gatk_CNNScoreVariants/output
    # FilterVariantTranches VCF files with FILTER tags
    o_gatk_FilterVariantTranches:
        type: File[]
        outputSource: gatk_FilterVariantTranches/output
    # CNN filtered VCF files - PASS only
    o_bcftools_view_filter_cnn:
        type: File[]
        outputSource: bcftools_view_filter_cnn/output
    o_bcftools_norm_cnn:
        type: File[]
        outputSource: bcftools_norm_cnn/output
    ### ANNOVAR outputs ###
    # CNN
    o_table_annovar_cnn_filtered_multianno_vcf:
        type: File[]
        outputSource: table_annovar_cnn_filtered/multianno_vcf
    o_table_annovar_cnn_filtered_multianno_txt:
        type: File[]
        outputSource: table_annovar_cnn_filtered/multianno_txt
    o_table_annovar_cnn_filtered_avinput:
        type: File[]
        outputSource: table_annovar_cnn_filtered/avinput
    # Hard filtering
    o_table_annovar_hard_filtered_multianno_vcf:
        type: File[]
        outputSource: table_annovar_hard_filtered/multianno_vcf
    o_table_annovar_hard_filtered_multianno_txt:
        type: File[]
        outputSource: table_annovar_hard_filtered/multianno_txt
    o_table_annovar_hard_filtered_avinput:
        type: File[]
        outputSource: table_annovar_hard_filtered/avinput

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
    ### Copy FASTQC output files - ZIP ###
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
    ### BWA MEM - Single end ###
    bwa_mem_single:
        run: ../wrappers/bwa-mem.cwl
        # scatterMethod: dotproduct
        scatter:
        - trimmed_fq_read1
        in:
            file_split: input_file_split
            sec_shorter_split_hits: bwa_mem_sec_shorter_split_hits
            num_threads: bwa_mem_num_threads
            ref_genome: reference_genome
            # R: 
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
    ### BWA MEM - Paired end ###
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
            # R: 
            trimmed_fq_read1: split_paired_read1_read2/reads_1
            trimmed_fq_read2: split_paired_read1_read2/reads_2
        out: [output]
    ### Collect and output all SAM files & generate filenames for next steps ###
    gather_bwa_sam_files:
        run: ../wrappers/gather-bwa-sam-files_single_v2.cwl
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
            isbam:
                valueFrom: $( true )
            collapsecigar: samtools_view_collapsecigar
            readsingroup: samtools_view_readsingroup
            uncompressed: samtools_view_uncompressed
            readtagtostrip: samtools_view_readtagtostrip
            input: gather_bwa_sam_files/total_sam_files
            readsquality: samtools_view_readsquality
            readswithbits: samtools_view_readswithbits
            # filter out using SAM FLAG values:
            readswithoutbits: samtools_view_readswithoutbits
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
                valueFrom: $(inputs.input.basename.split(".sam")[0].concat(".bam"))
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
                valueFrom: $(inputs.input_file.basename.split(".name.sorted.bam")[0].concat(".fixed.bam"))
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
                valueFrom: $(inputs.input.basename.split(".fixed.bam")[0].concat(".fixed.sorted.bam"))
        out: [sorted]
    ### Picard AddOrReplaceReadGroups - Assigns all the reads in a file to a single new read-group ###
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
                valueFrom: $(inputs.INPUT.basename.split(".fixed.sorted.bam")[0].concat(".fixed.sorted.uniq.rg.bam"))
            rgid: gather_bwa_sam_files/names_rgids
            rglb: gather_bwa_sam_files/names_rgplbs 
            rgpl: picard_addorreplacereadgroups_rgpl 
            rgpu: gather_bwa_sam_files/names_rgpus
            rgsm: 
                valueFrom: $(inputs.INPUT.basename.split(".fixed.sorted.bam")[0])
        out: [output]
    ### Picard MarkDuplicates ###
    picard_markduplicates:
        run: ../wrappers/picard-MarkDuplicates.cwl
        scatter: 
        - INPUT
        in:
            INPUT: picard_addorreplacereadgroups/output
            OUTPUT: 
                valueFrom: $(inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.bam")[0].concat(".fixed.sorted.uniq.rg.md.bam"))
            METRICS_FILE: 
                valueFrom: $(inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.bam")[0].concat(".md_metrics.txt"))
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
                valueFrom: $(inputs.input.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".align.stats.txt"))
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
            readswithbits: samtools_view_readswithbits
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
                valueFrom: $(inputs.input.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".count.total.txt"))
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
    # GATK BQSR
    ### samtools index ###
    samtools_index:
        run: ../wrappers/samtools-index.cwl
        scatter:
        - alignments
        in:
            alignments: picard_markduplicates/output
        out: [alignments_with_index]
    ### GATK BaseRecalibrator - GatherBQSRReports subworkflow ###
    gatk_bqsr_subworkflow:
        run: ../workflows/gatk-bqsr-subworkflow.cwl
        scatter: 
        - bqsr_INPUT
        in: 
            bqsr_INPUT: samtools_index/alignments_with_index
            bqsr_OUTPUT: 
                valueFrom: $(inputs.bqsr_INPUT.basename.split(".fixed.sorted.uniq.rg.md.bam")[0])
            bqsr_reference: reference_genome
            bqsr_known_sites_1: sub_bqsr_known_sites_1
            bqsr_known_sites_2: sub_bqsr_known_sites_2
            bqsr_known_sites_3: sub_bqsr_known_sites_3
            bqsr_intervals: gatk_splitintervals/intervalfiles
            bqsr_interval_padding: sub_bqsr_interval_padding
        out: 
        - o_gatk_gatherbqsrreports
    ### GATK ApplyBQSR ###
    gatk_applybqsr:
        run: ../wrappers/gatk-ApplyBQSR.cwl
        scatterMethod: dotproduct
        scatter: 
        - INPUT
        - bqsr_recal_file
        in: 
            INPUT: samtools_index/alignments_with_index 
            bqsr_recal_file: gatk_bqsr_subworkflow/o_gatk_gatherbqsrreports 
            OUTPUT: 
                valueFrom: $(inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".bqsr.bam"))
            reference: reference_genome
        out: [output]
    ### samtools index ###
    samtools_index_2:
        run: ../wrappers/samtools-index.cwl
        scatter:
        - alignments
        in:
            alignments: gatk_applybqsr/output
        out: [alignments_with_index]
    ### GATK HaplotypeCaller - Single-sample mode - without BQSR subworkflow ###
    gatk_haplotypecaller_subworkflow:
        run: ../workflows/gatk-haplotypecaller-subworkflow.cwl
        scatter: 
        - hf_INPUT
        in: 
            hf_INPUT: samtools_index/alignments_with_index
            hf_OUTPUT: 
                valueFrom: $(inputs.hf_INPUT.basename.split(".fixed.sorted.uniq.rg.md.bam")[0])
            hf_reference: reference_genome
            hf_intervals: gatk_splitintervals/intervalfiles
            hf_hc_native_pairHMM_threads: sub_hc_native_pairHMM_threads
            hf_hc_java_options: sub_hc_java_options
        out: 
        - o_gatk_MergeVCFs
    # GATK Hard-filtering
    ### GATK SelectVariants ###
    # SNPs
    gatk_SelectVariants_snps:
        run: ../wrappers/gatk-SelectVariants.cwl
        scatter: 
        - variant
        in:
            reference: reference_genome
            variant: gatk_haplotypecaller_subworkflow/o_gatk_MergeVCFs
            OUTPUT: 
              valueFrom: $(inputs.variant.basename.split(".vcf.gz")[0].concat(".snp.vcf"))
            variant_type: 
              valueFrom: $( "SNP" )
        out: [output]
    # INDELs
    gatk_SelectVariants_indels:
        run: ../wrappers/gatk-SelectVariants.cwl
        scatter: 
        - variant
        in:
            reference: reference_genome
            variant: gatk_haplotypecaller_subworkflow/o_gatk_MergeVCFs
            OUTPUT: 
              valueFrom: $(inputs.variant.basename.split(".vcf.gz")[0].concat(".indel.vcf"))
            variant_type: 
              valueFrom: $( "INDEL" )
        out: [output]
    ### GATK VariantFiltration ###
    # SNPs
    gatk_VariantFiltration_snps:
        run: ../wrappers/gatk-VariantFiltration.cwl
        scatter: 
        - variant
        in:
            reference: reference_genome
            variant: gatk_SelectVariants_snps/output
            window: VariantFiltration_window
            cluster: VariantFiltration_cluster
            filter_name: VariantFiltration_filter_name_snp
            filter: VariantFiltration_filter_snp
            OUTPUT: 
                valueFrom: $(inputs.variant.basename.split(".snp.vcf")[0].concat(".filtered.snp.vcf"))
        out: [output]
    # INDELs
    gatk_VariantFiltration_indels:
        run: ../wrappers/gatk-VariantFiltration.cwl
        scatter: 
        - variant
        in:
            reference: reference_genome
            variant: gatk_SelectVariants_indels/output
            window: VariantFiltration_window
            cluster: VariantFiltration_cluster
            filter_name: VariantFiltration_filter_name_indel
            filter: VariantFiltration_filter_indel
            OUTPUT: 
                valueFrom: $(inputs.variant.basename.split(".indel.vcf")[0].concat(".filtered.indel.vcf"))
        out: [output]
    ### bgzip - tabix ###
    # SNPs
    bgzip_snps:
        run: ../wrappers/bgzip.cwl
        scatter: 
        - input
        in: 
            input: gatk_VariantFiltration_snps/output
        out: [output]
    tabix_snps:
        run: ../wrappers/tabix.cwl
        scatter: 
        - input
        in: 
            input: bgzip_snps/output
        out: [output]
    # INDELs
    bgzip_indels:
        run: ../wrappers/bgzip.cwl
        scatter: 
        - input
        in: 
            input: gatk_VariantFiltration_snps/output
        out: [output]
    tabix_indels:
        run: ../wrappers/tabix.cwl
        scatter: 
        - input
        in: 
            input: bgzip_indels/output
        out: [output]
    ### bcftools concat - concatenate SNP and INDEL VCF files ###
    bcftools_concat:
        run: ../wrappers/bcftools-concat.cwl
        scatterMethod: dotproduct
        scatter: 
        - input1
        - input2
        in: 
            input1: tabix_snps/output
            input2: tabix_indels/output
            threads: bcftools_view_threads
            output_name: 
                valueFrom: $(inputs.input1.basename.split(".filtered.snp.vcf")[0].concat(".concat.vcf"))
        out: [output]
    ### bcftools view - filter variants ###
    bcftools_view_hard_filter:
        run: ../wrappers/bcftools-view.cwl
        scatter: 
        - input
        in: 
            input: bcftools_concat/output
            threads: bcftools_view_threads
            include: bcftools_view_include_hard_filters
            output_name: 
                valueFrom: $(inputs.input.basename.split(".concat.vcf")[0].concat(".bcftools.hard.filtered.vcf"))
        out: [output]
    ### bcftools norm ###
    bcftools_norm_hard_filter:
        run: ../wrappers/bcftools-norm.cwl
        scatter: 
        - input
        in: 
            input: bcftools_view_hard_filter/output
            threads: bcftools_norm_threads
            reference: reference_genome
            multiallelics: bcftoomls_norm_multiallelics
            output_type: 
                valueFrom: $( "v" )
        out: [output]
    ### ANNOVAR - variant annotation - Hard filtering ###
    table_annovar_hard_filtered:
        run: ../wrappers/table-annovar.cwl
        scatter: 
        - query_file
        in:
            query_file: bcftools_norm_hard_filter/output
            database_location: table_annovar_database_location
            build_over: table_annovar_build_over
            remove: table_annovar_remove
            protocol: table_annovar_protocol
            operation: table_annovar_operation
            na_string: table_annovar_na_string
            vcfinput: table_annovar_vcfinput
            otherinfo: table_annovar_otherinfo
            convert_arg: table_annovar_convert_arg
        out: [multianno_vcf, multianno_txt, avinput]

    # GATK CNN-filtering
    ### GATK CNNScoreVariants ###
    gatk_CNNScoreVariants:
        run: ../wrappers/gatk-CNNScoreVariants.cwl
        scatterMethod: dotproduct
        scatter:
        - variant
        - aligned_reads
        in: 
            reference: reference_genome
            variant: gatk_haplotypecaller_subworkflow/o_gatk_MergeVCFs
            aligned_reads: samtools_index_2/alignments_with_index
            OUTPUT: 
                valueFrom: $(inputs.variant.basename.split(".vcf.gz")[0].concat(".cnn.vcf"))
            tensor_type:
                valueFrom: $( "read_tensor" )
        out: [output]
    ### GATK FilterVariantTranches ###
    gatk_FilterVariantTranches:
        run: ../wrappers/gatk-FilterVariantTranches.cwl
        scatter:
        - variant
        in: 
            variant: gatk_CNNScoreVariants/output
            OUTPUT: 
                valueFrom: $(inputs.variant.basename.split(".cnn.vcf")[0].concat(".cnn_filtered.vcf"))
            resource_1: FilterVariantTranches_resource_1
            resource_2: FilterVariantTranches_resource_2
            resource_3: FilterVariantTranches_resource_3
        out: [output]
    ### bcftools view - filter variants ###
    bcftools_view_filter_cnn:
        run: ../wrappers/bcftools-view.cwl
        scatter: 
        - input
        in: 
            input: gatk_FilterVariantTranches/output
            threads: bcftools_view_threads
            include: bcftools_view_include_CNN_filters
            output_name: 
                valueFrom: $(inputs.input.basename.split(".cnn_filtered.vcf")[0].concat(".bcftools_cnn_filtered.vcf"))
        out: [output]
    ### bcftools norm ###
    bcftools_norm_cnn:
        run: ../wrappers/bcftools-norm.cwl
        scatter: 
        - input
        in: 
            input: bcftools_view_filter_cnn/output
            threads: bcftools_norm_threads
            reference: reference_genome
            multiallelics: bcftoomls_norm_multiallelics
            output_type: 
                valueFrom: $( "v" )
        out: [output]
    ### ANNOVAR - variant annotation - CNN ###
    table_annovar_cnn_filtered:
        run: ../wrappers/table-annovar.cwl
        scatter: 
        - query_file
        in:
            query_file: bcftools_norm_cnn/output
            database_location: table_annovar_database_location
            build_over: table_annovar_build_over
            remove: table_annovar_remove
            protocol: table_annovar_protocol
            operation: table_annovar_operation
            na_string: table_annovar_na_string
            vcfinput: table_annovar_vcfinput
            otherinfo: table_annovar_otherinfo
            convert_arg: table_annovar_convert_arg
        out: [multianno_vcf, multianno_txt, avinput]
