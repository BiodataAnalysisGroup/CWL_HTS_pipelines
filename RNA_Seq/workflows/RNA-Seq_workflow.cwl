class: Workflow
cwlVersion: v1.0

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

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
        default: false
    input_trimming_check:
        type: boolean?
        default: false
    premapping_input_check: 
        type: string
    ### Trim galore ###
    tg_quality: int
    tg_length: int
    tg_compression: boolean
    tg_do_not_compress: boolean
    tg_trim_suffix: string
    tg_strigency: int
    ### FASTA/Q Trimmer ###
    fastx_first_base_to_keep:
        type: int?
    fastx_last_base_to_keep:
        type: int?
    ### HISAT2 - mapping ###
    hisat2_num_of_threads:
        type: int
    hisat2_alignments_tailored_trans_assemb:
        type: boolean
    hisat2_idx_directory:
        type: Directory
    hisat2_idx_basename:
        type: string
    hisat2_known_splicesite_infile:
        type: File?
    ### samtools view ###
    samtools_view_isbam:
        type: boolean
        default: true
    samtools_view_collapsecigar:
        type: boolean
        default: false
    samtools_view_uncompressed:
        type: boolean
        default: false
    samtools_view_fastcompression:
        type: boolean
        default: false
    samtools_view_samheader:
        type: boolean
        default: false
    samtools_view_count:
        type: boolean
        default: false
    samtools_view_readswithoutbits: 
        type: int?
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
    ### samtools sort ###
    samtools_sort_compression_level:
        type: int?
    samtools_sort_threads:
        type: int?
    samtools_sort_memory:
        type: string?
    samtools_sort_sort_by_name:
        type: boolean?
        default: false
    ### Stringtie assembly of RNA-Seq alignments to transcripts ###
    stringtie_guide_gff:
        type: File
    stringtie_transcript_merge_mode:
        type: boolean
    stringtie_out_gtf:
        type: string
    stringtie_expression_estimation_mode:
        type: boolean
    stringtie_ballgown_table_files:
        type: boolean
    stringtie_cpus:
        type: int?
    stringtie_verbose:
        type: boolean?
        default: true
    stringtie_min_isoform_abundance:
        type: float?
    stringtie_junction_coverage:
        type: float?
    stringtie_min_read_coverage:
        type: float?
    stringtie_conservative_mode:
        default: false
        type: boolean?
    ### Ballgown - differential transcript expression analysis ###
    bg_phenotype_file:
        type: File
    bg_phenotype: 
        type: string
    bg_samples:
        type: string
    bg_timecourse:
        type: boolean?
    bg_feature:
        type: string?
    bg_measure:
        type: string?
    bg_confounders:
        type: string?
    bg_custom_model:
        type: boolean?
    bg_mod: 
        type: string?
    bg_mod0:
        type: string?
    ### featureCounts - counting reads ###
    featureCounts_number_of_threads:
        type: int?
        default: 16
    featureCounts_annotation_file:
        type: File
    featureCounts_output_file:
        type: string
    featureCounts_read_meta_feature_overlap:
        type: boolean?
    ### DESeq2 - differential gene expression analysis ###
    deseq2_metadata: 
        type: File
    deseq2_design:
        type: string
    deseq2_samples:
        type: string
    deseq2_min_sum_of_reads:
        type: int?
    deseq2_reference_level:
        type: string?
    deseq2_phenotype:
        type: string?
    deseq2_contrast: 
        type: boolean?
    deseq2_numerator: 
        type: string?
    deseq2_denominator:
        type: string?
    deseq2_lfcThreshold:
        type: float?
    deseq2_pAdjustMethod:
        type: string?
    deseq2_alpha:
        type: float?
    deseq2_parallelization:
        type: boolean?
    deseq2_cores:
        type: int?
    deseq2_transformation:
        type: string?
    deseq2_blind:
        type: boolean?
    deseq2_hypothesis:
        type: string?
    deseq2_reduced:
        type: string?
    deseq2_hidden_batch_effects:
        type: boolean?
    deseq2_hidden_batch_row_means:
        type: int?
    deseq2_hidden_batch_method:
        type: string?
    deseq2_variables:
        type: int?

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
    ### FASTA/Q Trimmer outputs ###
    o_fastx_trimmer_single:
        type: File[]
        outputSource: fastx_trimmer_single/output
    o_fastx_trimmer_paired:
        type: File[]
        outputSource: fastx_trimmer_paired/output
    ### HISAT2 outputs ###
    o_hisat2_for_single_reads_reports:
        type: File[]
        outputSource: hisat2_for_single_reads/output_stderr
    o_hisat2_for_paired_reads_reports:
        type: File[]
        outputSource: hisat2_for_paired_reads/output_stderr
    o_collect_hisat2_sam_files:
        type: File[]
        outputSource: collect_hisat2_sam_files/total_sam_files
    ### Samtools outputs ###
    o_samtools_view:
        type: File[]
        outputSource: samtools_view/output  
    o_samtools_sort:
        type: File[]
        outputSource: samtools_sort/sorted
    ### Stringtie outputs ###
    o_stringtie_transcript_assembly_gtf:
        type: File[]
        outputSource: stringtie_transcript_assembly/output_gtf
    o_stringtie_merge:
        type: File
        outputSource: stringtie_merge/output_gtf
    o_stringtie_expression_gtf:
        type: File[]
        outputSource: stringtie_expression/output_gtf
    o_stringtie_expression_outdir:
        type: Directory[]
        outputSource: stringtie_expression/outdir
    ### Ballgown outputs ###
    o_ballgown_de_results:
        type: File
        outputSource: ballgown_de/ballgown_de_results
    o_ballgown_object:
        type: File
        outputSource: ballgown_de/ballgown_object
    o_ballgown_de_custom_model:
        type: File?
        outputSource: ballgown_de/ballgown_de_custom_model
    ### featureCounts outputs ###
    o_featureCounts:
        type: File
        outputSource: featureCounts/output
    ### DESeq2 outputs ###
    o_deseq2_de_results:
        type: File
        outputSource: DESeq2_analysis/deseq2_de_results
    o_deseq2_dds_object:
        type: File
        outputSource: DESeq2_analysis/deseq2_dds_object
    o_deseq2_res_lfcShrink_object:
        type: File
        outputSource: DESeq2_analysis/deseq2_res_lfcShrink_object
    o_deseq2_transformed_object:
        type: File?
        outputSource: DESeq2_analysis/deseq2_transformed_object

steps:
    ### Retrieve FASTQ files from target directory (ExpressionTool) ###
    get_raw_files:
        run: ../wrappers/get-raw-files.cwl
        in:
            DIRECTORY: raw_files_directory
        out: [raw_files]
    ### Separate files - Generate file names (ExpressionTool) ###
    split_single_paired:
        run: ../wrappers/split-single-paired.cwl
        in:
            input_files: get_raw_files/raw_files
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            file_split_rev: input_file_split_rev
            qc_check: input_qc_check
            trimming_check: input_trimming_check
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
            - fastx_command_single
            - fastx_command_paired
    ### Trim_galore - FASTQ file trimming ###
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
    ### FASTQC - FASTQ file quality control ###
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
    ### FASTA/Q Trimmer - optional read processing ###
    fastx_trimmer_single:
        run: ../wrappers/fastx-trimmer.cwl
        scatter:
        - input_file
        in: 
            command: split_single_paired/fastx_command_single
            first_base_to_keep: fastx_first_base_to_keep
            last_base_to_keep: fastx_last_base_to_keep
            input_file: trim_galore_single/trim_galore 
        out: 
            - output
    fastx_trimmer_paired:
        run: ../wrappers/fastx-trimmer.cwl
        scatter:
        - input_file
        in: 
            command: split_single_paired/fastx_command_paired
            first_base_to_keep: fastx_first_base_to_keep
            last_base_to_keep: fastx_last_base_to_keep
            input_file: trim_galore_paired/trim_galore
        out: 
            - output
    ### Check inputs from selected steps and output file names for the next steps (ExpressionTool) ###
    check_for_fastx_and_produce_names:
        run: ../wrappers/check-selected-steps.cwl
        in: 
            input_check: premapping_input_check
            single_files: split_single_paired/single_files
            paired_files: split_single_paired/paired_files
            trim_galore_single: trim_galore_single/trim_galore
            trim_galore_paired: trim_galore_paired/trim_galore
            fastx_trimmer_single: fastx_trimmer_single/output
            fastx_trimmer_paired: fastx_trimmer_paired/output
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            file_split_rev: input_file_split_rev
        out: 
            - single_trim
            - single_hisat2_sam
            - paired_trim_1
            - paired_trim_2 
            - paired_hisat2_sam
    ### HISAT2 - Mapping to reference genome ###
    hisat2_for_single_reads:
        run: ../wrappers/hisat2.cwl
        scatterMethod: dotproduct
        scatter:
        - files_with_unpaired_reads
        - SAM_output
        in:
            num_of_threads: hisat2_num_of_threads
            alignments_tailored_trans_assemb: hisat2_alignments_tailored_trans_assemb
            idx_directory: hisat2_idx_directory
            idx_basename: hisat2_idx_basename
            files_with_unpaired_reads: check_for_fastx_and_produce_names/single_trim
            SAM_output: check_for_fastx_and_produce_names/single_hisat2_sam
        out: [output, output_stderr]
    hisat2_for_paired_reads:
        run: ../wrappers/hisat2.cwl
        scatterMethod: dotproduct
        scatter:
        - files_with_first_mates
        - files_with_second_mates
        - SAM_output
        in:
            num_of_threads: hisat2_num_of_threads
            known_splicesite_infile: hisat2_known_splicesite_infile
            alignments_tailored_trans_assemb: hisat2_alignments_tailored_trans_assemb
            idx_directory: hisat2_idx_directory
            idx_basename: hisat2_idx_basename
            files_with_first_mates: check_for_fastx_and_produce_names/paired_trim_1
            files_with_second_mates: check_for_fastx_and_produce_names/paired_trim_2
            SAM_output: check_for_fastx_and_produce_names/paired_hisat2_sam
        out: [output, output_stderr]
    ### Collect HISAT2 SAM files & create output file names for the remaining workflow steps (ExpressionTool) ###
    collect_hisat2_sam_files:
        run: ../wrappers/collect-hisat2-sam-files.cwl
        in:
            single_files: hisat2_for_single_reads/output
            paired_files: hisat2_for_paired_reads/output
        out: [total_sam_files]
    ### SAMtools ### 
    # SAM-to-BAM convertion
    samtools_view:
        run: ../wrappers/samtools-view.cwl
        scatter:
        - input
        in:
            input: collect_hisat2_sam_files/total_sam_files
            output_name:
                valueFrom: $( inputs.input.basename.split(".sam")[0].concat("_unsorted.bam") )
            threads: samtools_view_threads
            isbam: samtools_view_isbam
            collapsecigar: samtools_view_collapsecigar
            uncompressed: samtools_view_uncompressed
            fastcompression: samtools_view_fastcompression
            samheader: samtools_view_samheader
            count: samtools_view_count
            readswithoutbits: samtools_view_readswithoutbits
            readsingroup: samtools_view_readsingroup
            readtagtostrip: samtools_view_readtagtostrip
            readsquality: samtools_view_readsquality
            readswithbits: samtools_view_readswithbits
            cigar: samtools_view_cigar
            iscram: samtools_view_iscram
            randomseed: samtools_view_randomseed
            region: samtools_view_region
            readsinlibrary: samtools_view_readsinlibrary
        out: [output]
    # Sorting based on chromosomal coordinates
    samtools_sort:
        run: ../wrappers/samtools-sort.cwl
        scatter:
        - input
        in:
            input: samtools_view/output
            compression_level: samtools_sort_compression_level
            threads: samtools_sort_threads
            memory: samtools_sort_memory
            output_name: 
                valueFrom: $( inputs.input.basename.split("_unsorted.bam")[0].concat("_sorted.bam") )
            sort_by_name: samtools_sort_sort_by_name
        out: [sorted]
    ### Stringtie ### 
    # Transcript assembly
    stringtie_transcript_assembly:
        run: ../wrappers/stringtie.cwl
        scatter:
        - input_bam
        in:
            input_bam: samtools_sort/sorted
            guide_gff: stringtie_guide_gff
            sample_label: 
                valueFrom: $( inputs.input_bam.basename.split("_sorted.bam")[0] )
            cpus: stringtie_cpus
            out_gtf: 
                valueFrom: $( inputs.input_bam.basename.split("_sorted.bam")[0].concat(".gtf") )
            verbose: stringtie_verbose
            min_isoform_abundance: stringtie_min_isoform_abundance
            junction_coverage: stringtie_junction_coverage
            min_read_coverage: stringtie_min_read_coverage
            conservative_mode: stringtie_conservative_mode
        out: 
            - output_gtf    
    # Merge assembled transcripts
    stringtie_merge:
        run: ../wrappers/stringtie.cwl
        scatterMethod: dotproduct
        in:
            input_gtfs: stringtie_transcript_assembly/output_gtf
            transcript_merge_mode: stringtie_transcript_merge_mode
            guide_gff: stringtie_guide_gff
            cpus: stringtie_cpus
            verbose: stringtie_verbose
            out_gtf: stringtie_out_gtf
        out: 
            - output_gtf 
    # Generate transcript expression values
    stringtie_expression:
        run: ../wrappers/stringtie-for-ballgown.cwl
        scatter:
        - input_bam
        in:
            input_bam: samtools_sort/sorted
            guide_gff: stringtie_merge/output_gtf
            sample_label: 
                valueFrom: $( inputs.input_bam.basename.split("_sorted.bam")[0] )
            cpus: stringtie_cpus
            out_gtf: 
                valueFrom: $( inputs.input_bam.basename.split("_sorted.bam")[0].concat("_exprs.gtf") )
            verbose: stringtie_verbose
            min_isoform_abundance: stringtie_min_isoform_abundance
            junction_coverage: stringtie_junction_coverage
            min_read_coverage: stringtie_min_read_coverage
            conservative_mode: stringtie_conservative_mode
            outputdir: 
                valueFrom: $( inputs.input_bam.basename.split("_sorted.bam")[0] )
            expression_estimation_mode: stringtie_expression_estimation_mode
            ballgown_table_files: stringtie_ballgown_table_files
        out: 
            - output_gtf
            - outdir
    ### Ballgown - Differential transcript expression ###
    ballgown_de:
        run: ../wrappers/ballgown.cwl
        in: 
            stringtie_dirs: stringtie_expression/outdir
            phenotype_file: bg_phenotype_file
            phenotype: bg_phenotype
            samples: bg_samples
            timecourse: bg_timecourse
            feature: bg_feature
            measure: bg_measure
            confounders: bg_confounders
            custom_model: bg_custom_model
            mod: bg_mod
            mod0: bg_mod0
        out:
            - ballgown_de_results
            - ballgown_object
            - ballgown_de_custom_model
    ### featureCounts - Generation of read counts table ###
    featureCounts:
        run: ../wrappers/featureCounts.cwl
        in: 
            number_of_threads: featureCounts_number_of_threads
            annotation_file: featureCounts_annotation_file
            output_file: featureCounts_output_file
            inputFiles: samtools_sort/sorted
            read_meta_feature_overlap: featureCounts_read_meta_feature_overlap
        out: [output]
    ### DESeq2 - Differential gene expression ###
    DESeq2_analysis:
        run: ../wrappers/DESeq2.cwl
        in: 
            count_matrix: featureCounts/output
            metadata: deseq2_metadata
            samples: deseq2_samples
            design: deseq2_design
            min_sum_of_reads: deseq2_min_sum_of_reads
            reference_level: deseq2_reference_level
            phenotype: deseq2_phenotype
            contrast: deseq2_contrast
            numerator: deseq2_numerator
            denominator: deseq2_denominator
            lfcThreshold: deseq2_lfcThreshold
            pAdjustMethod: deseq2_pAdjustMethod
            alpha: deseq2_alpha
            parallelization: deseq2_parallelization
            cores: deseq2_cores
            transformation: deseq2_transformation
            blind: deseq2_blind
            hypothesis: deseq2_hypothesis
            reduced: deseq2_reduced
            hidden_batch_effects: deseq2_hidden_batch_effects
            hidden_batch_row_means: deseq2_hidden_batch_row_means
            hidden_batch_method: deseq2_hidden_batch_method
            variables: deseq2_variables
        out:
            - deseq2_de_results
            - deseq2_dds_object
            - deseq2_res_lfcShrink_object
            - deseq2_transformed_object
