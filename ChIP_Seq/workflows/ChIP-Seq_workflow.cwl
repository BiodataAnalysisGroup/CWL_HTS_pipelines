cwlVersion: v1.0
class: Workflow

requirements:
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

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
    ### Treatment & Control sample options ###
    input_treatment_samples:
        type: string[]
    input_control_samples:
        type: string[]?
    ### Trimmomatic ###
    # SE
    trimmomatic_se_threads: 
        type: int?
    trimmomatic_se_illuminaClip:
        type: string?
    trimmomatic_se_slidingWindow:
        type: string?
    trimmomatic_se_leading:
        type: int?
    trimmomatic_se_trailing:
        type: int?
    trimmomatic_se_minlen:
        type: int?
    # PE
    trimmomatic_pe_threads: 
        type: int?
    trimmomatic_pe_illuminaClip:
        type: string?
    trimmomatic_pe_slidingWindow:
        type: string?
    trimmomatic_pe_leading:
        type: int?
    trimmomatic_pe_trailing:
        type: int?
    trimmomatic_pe_minlen:
        type: int?
    ### HISAT2 - mapping ###
    hisat2_num_of_threads: 
        type: int?
    hisat2_idx_directory:
        type: Directory
    hisat2_idx_basename:
        type: string
    ### Samtools ###
    samtools_readswithoutbits: 
        type: int
        default: 4
    samtools_view_threads: 
        type: int?
    samtools_fixmate_threads:
        type: int?
    samtools_fixmate_output_format:
        type: string
        default: bam
    samtools_sort_compression_level:
        type: int?
    samtools_sort_threads:
        type: int
    samtools_sort_memory:
        type: string?
    samtools_markdup_threads: 
        type: int?    
    ### Deeptools ###
    # blacklist file
    blackListFile: 
        type: File?
    # multiBamSummary
    multiBamSummary_threads:
        type: int?
    # plotCorrelation
    plotCorrelation_numbers:
        type: boolean?
        default: true
    plotCorrelation_method:
        type: string?
    plotCorrelation_color:
        type: string?
    plotCorrelation_title:
        type: string?
    plotCorrelation_plotType:
        type: string?
    plotCorrelation_outFileName:
        type: string?
    # plotCoverage
    plotCoverage_threads:
        type: int?
    plotCoverage_plotFileFormat: 
        type: string?
    plotCoverage_outFileName: 
        type: string?
    # plotFingerprint
    plotFingerprint_plotFileFormat: 
        type: string?
    plotFingerprint_threads:
        type: int?
    plotFingerprint_outFileName: 
        type: string?
    # bamCoverage
    bamCoverage_normalizeUsing:
        type: string?
    bamCoverage_effective_genome_size: 
        type: long?
    bamCoverage_extendReads: 
        type: int?
    bamCoverage_threads:
        type: int?  
    # computeMatrix  
    computeMatrix_regions:
        type: File
    computeMatrix_threads: 
        type: int?
    computeMatrix_upstream:
        type: int?
    computeMatrix_downstream:
        type: int?
    computeMatrix_outputFile:
        type: string?
    computeMatrix_outFileSortedRegions:
        type: string?
    # plotHeatmap
    plotHeatmap_plotFileFormat: 
        type: string?
    plotHeatmap_outputFile: 
        type: string?
    ### MACS2 ###
    macs2_callpeak_bdg: 
        type: boolean?
    macs2_callpeak_gsize:
        type: string?
    macs2_callpeak_format:
        type: string?
    macs2_callpeak_broad:
        type: boolean?
    # --nomodel, --shift, --extsize
    macs2_callpeak_nomodel:
        type: boolean?
    macs2_shift: 
        type: int?
    macs2_extsize:  
        type: int?
    # p-/q-values
    macs2_pvalue: 
        type: float? 
    macs2_qvalue: 
        type: float?
    ### Differential binding ###
    ### ChIPQC & DiffBind ###
    metadata_table: 
        type: File
    ### ChIPQC ###
    ChIPQC_blacklist: 
        type: File?
    ChIPQC_annotation:
        type: string?
    ChIPQC_consensus:
        type: boolean?
    ChIPQC_bCount: 
        type: boolean? 
    ChIPQC_facetBy:
        type: string[]?
    ### DiffBind ###
    DiffBind_consensus: 
        type: string[]?
    DiffBind_minOverlap:
        type: 
        - int
        - float
        default: 2
    DiffBind_blacklist: 
        type: 
        - string
        - boolean
        - File
        default: true
    DiffBind_greylist: 
        type: 
        - string
        - boolean
        - File
        default: false
    DiffBind_cores:
        type: int?
    DiffBind_bParallel:
        type: boolean?
    DiffBind_normalization:
        type: string?
    DiffBind_library:
        type: string?
    DiffBind_background:
        type: boolean?
    DiffBind_design:
        type: string
        default: "~Condition"
    DiffBind_reorderMeta_factor:
        type: string[]?
    DiffBind_reorderMeta_value:
        type: string[]?
    DiffBind_retrieve_consensus: 
        type: boolean?
        default: true
    DiffBind_low_read_count_filter:
        type: int?
    DiffBind_filterFun:
        type: string?    
    ### ROSE ###
    rose_genome_build: 
        type: string
    rose_stitch_distance:
        type: int?
    rose_tss_distance:
        type: int?

outputs:
    ### Trimmomatic outputs ###
    o_trimmomatic_single_end_stderr:
        type: File[]
        outputSource: trimmomatic_single_end/stderr_log
    o_trimmomatic_single_end_fastq:
        type: File[]
        outputSource: trimmomatic_single_end/outFastq
    o_trimmomatic_paired_end_stderr:
        type: File[]
        outputSource: trimmomatic_paired_end/stderr_log
    o_trimmomatic_paired_end_fwd_paired:
        type: File[]
        outputSource: trimmomatic_paired_end/outFastq_fwd_paired
    o_trimmomatic_paired_end_fwd_unpaired:
        type: File[]
        outputSource: trimmomatic_paired_end/outFastq_fwd_unpaired
    o_trimmomatic_paired_end_rev_paired:
        type: File[]
        outputSource: trimmomatic_paired_end/outFastq_rev_paired
    o_trimmomatic_paired_end_rev_unpaired:
        type: File[]
        outputSource: trimmomatic_paired_end/outFastq_rev_unpaired
   ### FASTQC outputs ###
    # HTML
    o_fastqc_raw_html:
        type: File[]?
        outputSource: rename_fastqc_raw_html/renamed_file
    o_fastqc_single_html:
        type: File[]?
        outputSource: rename_fastqc_single_html/renamed_file
    o_fastqc_paired_html_fwd:
        type: File[]?
        outputSource: rename_fastqc_paired_html_fwd/renamed_file
    o_fastqc_paired_html_rev:
        type: File[]?
        outputSource: rename_fastqc_paired_html_rev/renamed_file
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
    ### HISAT2 outputs ###
    o_hisat2_for_single_reads_sam:
        type: File[]
        outputSource: hisat2_for_single_reads/output
    o_hisat2_for_single_reads_stderr:
        type: File[]
        outputSource: hisat2_for_single_reads/output_stderr
    o_hisat2_for_paired_reads_sam:
        type: File[]
        outputSource: hisat2_for_paired_reads/output
    o_hisat2_for_paired_reads_stderr:
        type: File[]
        outputSource: hisat2_for_paired_reads/output_stderr
    ### samtools outputs ###
    o_samtools_sort_by_name:
        type: File[]
        outputSource: samtools_sort_by_name/sorted
    # Add ms (mate score) tags
    o_samtools_fixmate:
        type: File[]
        outputSource: samtools_fixmate/output
    # Sort by coordinates
    o_samtools_sort:
        type: File[]
        outputSource: samtools_sort/sorted
    # Remove duplicates (samtools markdup)
    o_samtools_markdup:
        type: File[]
        outputSource: samtools_markdup/output
    # Index bam
    o_samtools_index:
        type: File[]
        outputSource: samtools_index/alignments_with_index
    ### Deeptools outputs ###
    # Multibam file
    o_multiBamSummary_file:
        type: File
        outputSource: multiBamSummary_file/outNpz
    # plotCorrelation image file
    o_plotCorrelation_file:
        type: File
        outputSource: plotCorrelation_file/outImage
    # plotCoverage image files
    o_plotCoverage_file:
        type: File
        outputSource: plotCoverage_file/outImage
    # plotFingerprint image files
    o_plotFingerprint_file:
        type: File
        outputSource: plotFingerprint_file/outImage
    # bamCoverage normalization
    o_bamCoverage_norm:
        type: File[]
        outputSource: bamCoverage_norm/bigwig
    o_computeMatrix_matrix:
        type: File
        outputSource: computeMatrix/outMatrix
    o_computeMatrix_regions:
        type: File
        outputSource: computeMatrix/outRegions
    o_plotHeatmap:
        type: File
        outputSource: plotHeatmap/outHeatmap
    ### MACS2 outputs ### 
    o_macs2_call_peaks_narrowPeak:
        type: File[]?
        outputSource: macs2_call_peaks/narrowPeak
    o_macs2_call_peaks_xls:
        type: File[]?
        outputSource: macs2_call_peaks/xls
    o_macs2_call_peaks_bed:
        type: File[]?
        outputSource: macs2_call_peaks/bed
    o_macs2_call_peaks_lambda:
        type: File[]?
        outputSource: macs2_call_peaks/lambda
    o_macs2_call_peaks_pileup:
        type: File[]?
        outputSource: macs2_call_peaks/pileup
    o_macs2_call_peaks_broadPeak:
        type: File[]?
        outputSource: macs2_call_peaks/broadPeak
    o_macs2_call_peaks_gappedPeak:
        type: File[]?
        outputSource: macs2_call_peaks/gappedPeak
    o_macs2_call_peaks_model_r:
        type: File[]?
        outputSource: macs2_call_peaks/model_r
    o_macs2_call_peaks_cutoff:
        type: File[]?
        outputSource: macs2_call_peaks/cutoff
    ### Consensus peaks and table with counts ###
    o_total_peaks_table:
        type: File
        outputSource: total_peaks_table/output
    o_sort_peaks_table:
        type: File
        outputSource: sort_peaks_table/output
    o_bedtools_merge:
        type: File
        outputSource: bedtools_merge/output
    o_bedtools_intersect:
        type: File
        outputSource: bedtools_merge/output
    o_exclude_black_list_regions:
        type: File
        outputSource: exclude_black_list_regions/output
    o_bedtools_coverage:
        type: File[]
        outputSource: bedtools_coverage/output
    o_printf_header_samples:
        type: File
        outputSource: printf_header_samples/output
    o_paste_content_1:
        type: File
        outputSource: paste_content_1/output
    o_paste_content_2:
        type: File
        outputSource: paste_content_2/output
    o_append_files:
        type: File
        outputSource: append_files/output
    ### ChIPQC MACS output ###
    o_ChIPQC_macs_ChIPQCexperiment:
        type: File
        outputSource: ChIPQC_macs/ChIPQCexperiment
    o_ChIPQC_macs_outdir:
        type: Directory
        outputSource: ChIPQC_macs/outdir
    o_ChIPQC_macs_ChIPQCreport:
        type: File?
        outputSource: ChIPQC_macs/ChIPQCreport
    ### DiffBind MACS output ###
    o_DiffBind_macs_diffbind_results:
        type: File
        outputSource: DiffBind_macs/diffbind_results
    o_DiffBind_macs_correlation_heatmap:
        type: File
        outputSource: DiffBind_macs/correlation_heatmap
    o_DiffBind_macs_diffbind_consensus:
        type: File?
        outputSource: DiffBind_macs/diffbind_consensus
    o_DiffBind_macs_diffbind_normalized_counts:
        type: File?
        outputSource: DiffBind_macs/diffbind_normalized_counts
    o_DiffBind_macs_diffbind_dba_object:
        type: File?
        outputSource: DiffBind_macs/diffbind_dba_object
    ### ROSE-related outputs ###
    # exclude_black_list_regions_narrowPeak
    o_exclude_black_list_regions_narrowPeak:
        type: File[]
        outputSource: exclude_black_list_regions_narrowPeak/output
    # bed to rose-compatible gff
    o_bed_to_rose_gff_conversion:
        type: File[]
        outputSource: bed_to_rose_gff_conversion/output 
    # ROSE Super-Enhancers output  
    o_rose_main_gff_dir_outputs:
        type: {
            "type": "array",
            "items": {"type" : "array", "items" : "File"}
        }
        outputSource: rose_main/gff_dir_outputs
    o_rose_main_mappedGFF_dir_outputs:
        type: {
            "type": "array",
            "items": {"type" : "array", "items" : "File"}
        }
        outputSource: rose_main/mappedGFF_dir_outputs
    o_rose_main_STITCHED_ENHANCER_REGION_MAP:
        type: File[]?
        outputSource: rose_main/STITCHED_ENHANCER_REGION_MAP
    o_rose_main_AllEnhancers_table:
        type: File[]?
        outputSource: rose_main/AllEnhancers_table
    o_rose_main_SuperEnhancers_table:
        type: File[]?
        outputSource: rose_main/SuperEnhancers_table
    o_rose_main_Plot_points:
        type: File[]?
        outputSource: rose_main/Plot_points
    o_rose_main_Enhancers_withSuper:
        type: File[]?
        outputSource: rose_main/Enhancers_withSuper
    # enhancer_bed_processing
    o_enhancer_bed_processing:
        type: File[]?
        outputSource: enhancer_bed_processing/output
    # ChIPQC ROSE output
    o_ChIPQC_rose_ChIPQCexperiment:
        type: File
        outputSource: ChIPQC_rose/ChIPQCexperiment
    o_ChIPQC_rose_outdir:
        type: Directory?
        outputSource: ChIPQC_rose/outdir
    o_ChIPQC_rose_ChIPQCreport:
        type: File
        outputSource: ChIPQC_rose/ChIPQCreport
    # DiffBind ROSE output
    o_DiffBind_rose_diffbind_results:
        type: File
        outputSource: DiffBind_rose/diffbind_results
    o_DiffBind_rose_correlation_heatmap:
        type: File
        outputSource: DiffBind_rose/correlation_heatmap
    o_DiffBind_rose_diffbind_consensus:
        type: File?
        outputSource: DiffBind_rose/diffbind_consensus
    o_DiffBind_rose_diffbind_normalized_counts:
        type: File?
        outputSource: DiffBind_rose/diffbind_normalized_counts
    o_DiffBind_rose_diffbind_dba_object:
        type: File?
        outputSource: DiffBind_rose/diffbind_dba_object

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
            - paired_files_fwd
            - paired_files_rev
            - fastqc_for_raw
            - fastqc_for_single
            - fastqc_for_paired
            - single_files_basenames
            - paired_files_basenames
            - single_files_sam
            - paired_files_sam
            - cp_command_raw
            - cp_command_single
            - cp_command_paired
            - trimmomatic_command_single
            - trimmomatic_command_paired
    ### Trimmomatic - trimming ###
    trimmomatic_single_end:
        run: ../wrappers/Trimmomatic_SE.cwl
        # scatterMethod: dotproduct
        scatter:
            - input_fastq
        in: 
            input_fastq: split_single_paired/single_files
            trimm_se_threads: trimmomatic_se_threads
            command: split_single_paired/trimmomatic_command_single
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            illuminaClip: trimmomatic_se_illuminaClip
            slidingWindow: trimmomatic_se_slidingWindow
            leading: trimmomatic_se_leading
            trailing: trimmomatic_se_trailing
            minlen: trimmomatic_se_minlen
        out:
            - stderr_log
            - outFastq
    trimmomatic_paired_end:
        run: ../wrappers/Trimmomatic_PE.cwl
        scatter:
            - input_fastq_fwd
            - input_fastq_rev
        scatterMethod: dotproduct
        in: 
            input_fastq_fwd: split_single_paired/paired_files_fwd
            input_fastq_rev: split_single_paired/paired_files_rev
            trimm_pe_threads: trimmomatic_pe_threads
            command: split_single_paired/trimmomatic_command_paired
            file_split: input_file_split
            file_split_fwd_single: input_file_split_fwd_single
            file_split_rev: input_file_split_rev
            illuminaClip: trimmomatic_pe_illuminaClip
            slidingWindow: trimmomatic_pe_slidingWindow
            leading: trimmomatic_pe_leading
            trailing: trimmomatic_pe_trailing
            minlen: trimmomatic_pe_minlen
        out:
            - stderr_log
            - outFastq_fwd_paired
            - outFastq_fwd_unpaired
            - outFastq_rev_paired
            - outFastq_rev_unpaired
    ### FASTQC - Quality control ###
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
            input_files: trimmomatic_single_end/outFastq
        out:
            - html_file
            - zipped_file
    fastqc_paired_trimmed_fwd:
        run: ../wrappers/fastqc.cwl
        in:
            command: split_single_paired/fastqc_for_paired
            input_files: trimmomatic_paired_end/outFastq_fwd_paired
        out:
            - html_file
            - zipped_file
    fastqc_paired_trimmed_rev:
        run: ../wrappers/fastqc.cwl
        in:
            command: split_single_paired/fastqc_for_paired
            input_files: trimmomatic_paired_end/outFastq_rev_paired
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
        run: ../wrappers/cp_paired.cwl
        in:
            command: split_single_paired/cp_command_paired
            input_files_fwd: fastqc_paired_trimmed_fwd/zipped_file
            input_files_rev: fastqc_paired_trimmed_rev/zipped_file
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
    rename_fastqc_paired_html_fwd:
        run: ../wrappers/rename.cwl
        in:
            input_file: fastqc_paired_trimmed_fwd/html_file
            run_type:
                valueFrom: $( "fastqc_paired_trimmed_" )
        out:
            - renamed_file
        scatter:
            - input_file   
    rename_fastqc_paired_html_rev:
        run: ../wrappers/rename.cwl
        in:
            input_file: fastqc_paired_trimmed_rev/html_file
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
            input_paired_fwd: split_single_paired/paired_files_fwd
            input_paired_rev: split_single_paired/paired_files_rev
            trimming_single: trimmomatic_single_end/outFastq
            trimming_paired_fwd: trimmomatic_paired_end/outFastq_fwd_paired
            trimming_paired_rev: trimmomatic_paired_end/outFastq_rev_paired
        out:
        - single_files
        - paired_files_fwd
        - paired_files_rev
    ### HISAT2 ###
    hisat2_for_single_reads:
        run: ../wrappers/hisat2_v2.cwl 
        scatterMethod: dotproduct
        scatter:
            - files_with_unpaired_reads
            - SAM_output
            - stderr_report
        in:
            num_of_threads: hisat2_num_of_threads
            idx_directory: hisat2_idx_directory
            idx_basename: hisat2_idx_basename
            files_with_unpaired_reads: check_trimming/single_files
            SAM_output: split_single_paired/single_files_sam
            stderr_report: split_single_paired/single_files_basenames
        out: 
            - output
            - output_stderr
    hisat2_for_paired_reads:
        run: ../wrappers/hisat2_v2.cwl
        scatterMethod: dotproduct
        scatter:
            - files_with_first_mates
            - files_with_second_mates
            - SAM_output
            - stderr_report
        in:
            num_of_threads: hisat2_num_of_threads
            idx_directory: hisat2_idx_directory
            idx_basename: hisat2_idx_basename
            files_with_first_mates: check_trimming/paired_files_fwd
            files_with_second_mates: check_trimming/paired_files_rev
            SAM_output: split_single_paired/paired_files_sam
            stderr_report: split_single_paired/paired_files_basenames
        out: 
            - output
            - output_stderr
    ### Collect HISAT2 SAM files ###
    collect_hisat2_sam_files:
        run: ../wrappers/collect-hisat2-sam-files.cwl
        in:
            single_files: hisat2_for_single_reads/output
            paired_files: hisat2_for_paired_reads/output           
        out:
            - total_sam_files
    ### Samtools view ###
    samtools_view:
        run: ../wrappers/samtools-view.cwl
        scatter: 
        - input
        in: 
            input: collect_hisat2_sam_files/total_sam_files 
            output_name:
                valueFrom: $( inputs.input.basename.split(".sam")[0].concat(".bam.tmp") )
            readswithoutbits: samtools_readswithoutbits
            samheader:  
                valueFrom: $( true )
            threads: samtools_view_threads
        out: [output]
    ### Samtools sort by read names ###
    samtools_sort_by_name:
        run: ../wrappers/samtools-sort.cwl
        scatter:
        - input
        in:
            compression_level: samtools_sort_compression_level 
            threads: samtools_sort_threads 
            memory: samtools_sort_memory 
            input: samtools_view/output
            output_name: 
                valueFrom: $( inputs.input.basename.split(".bam.tmp")[0].concat(".name.sorted.bam") )
            sort_by_name:
                valueFrom: $( true )
        out: [sorted]
    ### Samtools fixmate ###
    samtools_fixmate:
        run: ../wrappers/samtools-fixmate.cwl
        scatter:
        - input_file
        in:
            threads: samtools_fixmate_threads
            output_format: samtools_fixmate_output_format
            input_file: samtools_sort_by_name/sorted
            output_file_name: 
                valueFrom: $( inputs.input_file.basename.split(".name.sorted.bam")[0].concat("_unsorted.bam") )
        out: [output]
    ### sort BAM by CHR coordinates ###
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
                valueFrom: $( inputs.input.basename.split("_unsorted.bam")[0].concat("_sorted.bam") )
        out: [sorted]
    ### Remove duplicates ###
    samtools_markdup:
        run: ../wrappers/samtools-markdup.cwl
        scatter: 
            - input
        in: 
            threads: samtools_markdup_threads
            remove_duplicates: 
                valueFrom: $( true )
            input: samtools_sort/sorted
            output_name: 
                valueFrom: $( inputs.input.basename.split("_sorted.bam")[0].concat("_markdup.bam") )
        out: [output]
    ### Index BAM files ###
    samtools_index:
        run: ../wrappers/samtools-index.cwl
        scatter:
            - alignments
        in: 
            alignments: samtools_markdup/output
        out: [alignments_with_index]
    
    ### Deeptools analyses ###
    ### Create multibam file ###
    multiBamSummary_file:
        run: ../wrappers/multiBamSummary.cwl
        in: 
            bam: samtools_index/alignments_with_index
            threads: multiBamSummary_threads
            blackListFileName: blackListFile
        out: [outNpz]
    ### Plot heatmap ###
    plotCorrelation_file:
        run: ../wrappers/plotCorrelation.cwl
        in: 
            input: multiBamSummary_file/outNpz
            numbers: plotCorrelation_numbers
            method: plotCorrelation_method
            color: plotCorrelation_color
            title: plotCorrelation_title
            plotType: plotCorrelation_plotType
            outFileName: plotCorrelation_outFileName
        out: [outImage]
    ## Plot read coverage ###
    plotCoverage_file:
        run: ../wrappers/plotCoverage.cwl
        in:
            bam: samtools_index/alignments_with_index
            threads: plotCoverage_threads
            skipZeros: 
                valueFrom: $( true )
            plotFileFormat: plotCoverage_plotFileFormat
            outFileName: plotCoverage_outFileName
            blackListFileName: blackListFile
        out: [outImage]
    ### Plot fingerprint ###
    plotFingerprint_file:
        run: ../wrappers/plotFingerprint.cwl
        in:
            bam: samtools_index/alignments_with_index
            threads: plotFingerprint_threads
            plotFileFormat: plotFingerprint_plotFileFormat
            outFileName: plotFingerprint_outFileName
            blackListFileName: blackListFile
        out: [outImage]
    ### Convert bam to bigwig and normalize ###
    # Normalize using BPM (bins per million mapped reads) OR
    # Normalize using RPGC (reads per genome coverage, effective genome size MM10) OR
    # Normalize using RPKM (reads Per Kilobase per Million mapped reads)
    bamCoverage_norm:
        run: ../wrappers/bamCoverage.cwl
        scatter: 
            - bam
        in: 
            bam: samtools_index/alignments_with_index
            effective_genome_size: bamCoverage_effective_genome_size
            normalizeUsing: bamCoverage_normalizeUsing
            extendReads: bamCoverage_extendReads
            threads: bamCoverage_threads
            blackListFileName: blackListFile
        out: [bigwig]
    ### Computematrix ###
    computeMatrix:
        run: ../wrappers/computeMatrix.cwl
        in: 
            bw: bamCoverage_norm/bigwig
            regions: computeMatrix_regions
            threads: computeMatrix_threads
            upstream: computeMatrix_upstream
            downstream: computeMatrix_downstream
            outputFile: computeMatrix_outputFile
            outFileSortedRegions: computeMatrix_outFileSortedRegions
            blackListFileName: blackListFile
        out: [outMatrix, outRegions]
    ### Plot heatmap ###
    plotHeatmap:
        run: ../wrappers/plotHeatmap.cwl
        in:
            matrix: computeMatrix/outMatrix
            plotFileFormat: plotHeatmap_plotFileFormat
            outputFile: plotHeatmap_outputFile
        out: [outHeatmap]
    
    ### PEAK CALLING ###
    ### Separate control and treatment samples for MACS2 ###
    separate_control_treatment_files:
        run: ../wrappers/separate_control_treatment_files.cwl
        in: 
            treatment_samples: input_treatment_samples
            control_samples: input_control_samples
            aligned_files: samtools_index/alignments_with_index
        out:
        - treatment_files
        - control_files
    ### Call peaks ###
    macs2_call_peaks:
        run: ../wrappers/macs2-callpeak.cwl
        in:
            bdg: macs2_callpeak_bdg
            treatment: separate_control_treatment_files/treatment_files
            control: separate_control_treatment_files/control_files
            f: macs2_callpeak_format
            gsize: macs2_callpeak_gsize
            broad: macs2_callpeak_broad
            # --nomodel, --shift, --extsize
            nomodel: macs2_callpeak_nomodel
            shift: macs2_shift
            extsize: macs2_extsize
            # p-/q-values
            p: macs2_pvalue
            q: macs2_qvalue
        scatter:
        - treatment
        - control
        scatterMethod: dotproduct
        out:
        - narrowPeak
        - xls
        - bed
        - lambda
        - pileup
        - broadPeak
        - gappedPeak
        - model_r
        - cutoff
    
    ### Create consensus peaks and table with counts for custom analyses ###
    ### Create table with total peaks ###
    total_peaks_table:
        run: ../wrappers/table_total_peaks.cwl
        in: 
            input: macs2_call_peaks/narrowPeak
            consensus_bed:
                valueFrom: $( "narrowPeaks_consensus.bed" )
        out: [output]
    ### Sort peaks ###
    sort_peaks_table:
        run: ../wrappers/sort.cwl
        in: 
            input: total_peaks_table/output
            key1: 
                valueFrom: $("1,1")
            key2: 
                valueFrom: $("2,2n")
            output_name:
                valueFrom: $( "narrowPeaks_consensus_sorted.bed" )
        out: [output]
    ### Intersect peaks to consensus ###
    bedtools_merge:
        run: ../wrappers/bedtools-merge.cwl
        in: 
            input: sort_peaks_table/output
            output_name:
                valueFrom: $( "narrowPeaks_consensus_merged.bed" )
        out: [output]
    ### Exclude blacklisted regions ###
    exclude_black_list_regions:
        run: ../wrappers/bedtools-intersect.cwl
        in:
            feature_a: bedtools_merge/output
            feature_b: blackListFile
            no_overlap: 
                valueFrom: $( true )
            outputFile:
                valueFrom: $( "narrowPeaks_consensus_filtered.bed" )
        out: [output]
    ### Calculate read coverage for every peak on the consensus ###
    bedtools_coverage:
        run: ../wrappers/bedtools-coverage.cwl
        scatter: 
        - input_file
        in: 
            reference_file: exclude_black_list_regions/output
            input_file: separate_control_treatment_files/treatment_files
            count_of_overlaps: 
                valueFrom: $( true )
        out: [output]
    extract_counts: 
        run: ../wrappers/awk.cwl
        scatter: 
         - input_file
        in: 
            print: 
                valueFrom: '{print $4}'
            input_file: bedtools_coverage/output
            split_string: 
                valueFrom: $("_markdup_counts.bed")
            output_suffix: 
                valueFrom: $(".count")
        out: [output]
    extract_peaks:
        run: ../wrappers/awk.cwl
        in: 
            print: 
                valueFrom: '{print $1,$2,$3}'
            input_file: exclude_black_list_regions/output
            split_string: 
                valueFrom: $(".bed")
            output_suffix: 
                valueFrom: $(".coords")
        out: [output]
    ### Merge final counts with headers ###
    printf_header_samples:
        run: ../wrappers/printf.cwl
        in:
            input_files: bedtools_coverage/output
            output_name:
                valueFrom: $( "final_consensus_count_samples.bed" )
        out: [output]
    paste_content_1:
        run: ../wrappers/paste.cwl
        in: 
            input_files: extract_counts/output
            output_name:
                valueFrom: $( "consensus_count.txt" )
        out: [output] 
    paste_content_2:
        run: ../wrappers/paste.cwl
        in: 
            coordinates: extract_peaks/output
            counts: paste_content_1/output
            output_name:
                valueFrom: $( "consensus_count_with_coordinates.txt" )
        out: [output] 
    append_files:
        run: ../wrappers/awk.cwl
        in: 
            print: 
                valueFrom: '{print $0}'
            input_file: printf_header_samples/output
            second_file: paste_content_2/output
            split_string: 
                valueFrom: $("_samples.bed")
            output_suffix: 
                valueFrom: $(".txt")
        out: [output]

    ### QC and Differential Binding analysis for MACS2 peaks ###
    ### ChIPQC ###
    ChIPQC_macs:
        run: ../wrappers/ChIPQC.cwl
        in: 
            treatmentBAM: separate_control_treatment_files/treatment_files
            controlBAM: separate_control_treatment_files/control_files
            peaks: macs2_call_peaks/narrowPeak
            peakCaller: 
                valueFrom: $("narrow")
            annotation: ChIPQC_annotation
            metadata: metadata_table
            consensus: ChIPQC_consensus
            bCount: ChIPQC_bCount
            blacklist: ChIPQC_blacklist
            facetBy: ChIPQC_facetBy
            outputdir: 
                valueFrom: $("ChIPQC_macs_HTML_report")
            chipqc_experiment: 
                valueFrom: $("ChIPQC_macs_experiment.rda")
            chipqc_report: 
                valueFrom: $("ChIPQC_macs_report.csv")
        out:
        - ChIPQCexperiment
        - outdir
        - ChIPQCreport
    ### DiffBind ###
    DiffBind_macs: 
        run: ../wrappers/DiffBind.cwl
        in: 
            treatmentBAM: separate_control_treatment_files/treatment_files
            controlBAM: separate_control_treatment_files/control_files
            peaks: macs2_call_peaks/narrowPeak
            peakCaller: 
                valueFrom: $("narrow")
            metadata: metadata_table
            consensus: DiffBind_consensus
            minOverlap: DiffBind_minOverlap
            blacklist: DiffBind_blacklist
            greylist: DiffBind_greylist
            cores: DiffBind_cores
            bParallel: DiffBind_bParallel
            normalization: DiffBind_normalization
            library: DiffBind_library
            background: DiffBind_background
            design: DiffBind_design
            reorderMeta_factor: DiffBind_reorderMeta_factor
            reorderMeta_value: DiffBind_reorderMeta_value
            retrieve_consensus: DiffBind_retrieve_consensus
            low_read_count_filter: DiffBind_low_read_count_filter
            diffbind_filterFun: DiffBind_filterFun
            diffbind_results_filename: 
                valueFrom: $("DiffBind_macs_diff_peaks.tsv")
            correlation_heatmap_filename:
                valueFrom: $("DiffBind_macs_diff_peaks_corr_heatmap.tif")
            diffbind_consensus_filename:
                valueFrom: $("DiffBind_macs_diff_consensus_peaks.bed")
        out:
        - diffbind_results
        - correlation_heatmap
        - diffbind_consensus
        - diffbind_normalized_counts
        - diffbind_dba_object
    
    ### CALL SUPER ENHANCERS ###
    ### Exclude blacklisted regions for ROSE input ###
    exclude_black_list_regions_narrowPeak:
        run: ../wrappers/bedtools-intersect-narrowPeak.cwl
        scatter: 
        - feature_a
        in:
            feature_a: macs2_call_peaks/narrowPeak
            feature_b: blackListFile
            no_overlap: 
                valueFrom: $( true )
        out: [output]
    ### Convert to ROSE gff format ###
    bed_to_rose_gff_conversion:
        run: ../wrappers/bed_to_rose_gff.cwl
        scatter: 
        - input
        in: 
            input: exclude_black_list_regions_narrowPeak/output
        out: [output]
    ### Identify Super Enhancers ###
    rose_main:
        run: ../wrappers/rose.cwl
        scatterMethod: dotproduct
        scatter: 
            - macs_peaks
            - ranking_bam
            - control_bam
        in: 
            macs_peaks: bed_to_rose_gff_conversion/output
            ranking_bam: separate_control_treatment_files/treatment_files
            genome_build: rose_genome_build
            stitch_distance: rose_stitch_distance
            tss_distance: rose_tss_distance
            control_bam: separate_control_treatment_files/control_files
        out:
            - gff_dir_outputs
            - mappedGFF_dir_outputs
            - STITCHED_ENHANCER_REGION_MAP
            - AllEnhancers_table
            - SuperEnhancers_table
            - Enhancers_withSuper
            - Plot_points
    # QC and Differential Binding analysis for ROSE super-enhancers
    ### Adjust BED files with ROSE Enhancer and Super-Enhancer information ###
    enhancer_bed_processing:
        run: ../wrappers/enhancer_bed_processing.cwl
        scatter: 
        - input
        in: 
            input: rose_main/Enhancers_withSuper
        out: [output]
    ### ChIPQC ###
    ChIPQC_rose:
        run: ../wrappers/ChIPQC.cwl
        in: 
            treatmentBAM: separate_control_treatment_files/treatment_files
            controlBAM: separate_control_treatment_files/control_files
            peaks: enhancer_bed_processing/output
            peakCaller: 
                valueFrom: $("bed") 
            annotation: ChIPQC_annotation
            metadata: metadata_table
            consensus: ChIPQC_consensus
            bCount: ChIPQC_bCount
            blacklist: ChIPQC_blacklist
            facetBy: ChIPQC_facetBy
            outputdir: 
                valueFrom: $("ChIPQC_rose_HTML_report")
            chipqc_experiment: 
                valueFrom: $("ChIPQC_rose_experiment.rda")
            chipqc_report: 
                valueFrom: $("ChIPQC_rose_report.csv")
        out:
        - ChIPQCexperiment
        - outdir
        - ChIPQCreport
    ### DiffBind ###
    DiffBind_rose: 
        run: ../wrappers/DiffBind.cwl
        in: 
            treatmentBAM: separate_control_treatment_files/treatment_files
            controlBAM: separate_control_treatment_files/control_files
            peaks: enhancer_bed_processing/output
            peakCaller: 
                valueFrom: $("bed") 
            metadata: metadata_table
            consensus: DiffBind_consensus
            minOverlap: DiffBind_minOverlap
            blacklist: DiffBind_blacklist
            greylist: DiffBind_greylist
            cores: DiffBind_cores
            bParallel: DiffBind_bParallel
            normalization: DiffBind_normalization
            library: DiffBind_library
            background: DiffBind_background
            design: DiffBind_design
            reorderMeta_factor: DiffBind_reorderMeta_factor
            reorderMeta_value: DiffBind_reorderMeta_value
            retrieve_consensus: DiffBind_retrieve_consensus
            low_read_count_filter: DiffBind_low_read_count_filter
            diffbind_filterFun: DiffBind_filterFun
            diffbind_results_filename: 
                valueFrom: $("DiffBind_rose_diff_peaks.tsv")
            correlation_heatmap_filename:
                valueFrom: $("DiffBind_rose_diff_peaks_corr_heatmap.tif")
            diffbind_consensus_filename:
                valueFrom: $("DiffBind_rose_diff_consensus_peaks.bed")
        out:
        - diffbind_results
        - correlation_heatmap
        - diffbind_consensus
        - diffbind_normalized_counts
        - diffbind_dba_object
