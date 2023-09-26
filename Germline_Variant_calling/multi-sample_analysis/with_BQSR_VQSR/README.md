# CWL-based (multi-sample) workflow for germline variant calling

## Description

A CWL-based pipeline for calling small germline variants, namely SNPs and small INDELs, by processing data from Whole-genome Sequencing (WGS) or Targeted Sequencing (e.g., Whole-exome sequencing; WES) experiments. A pre-configured YAML template, based on validation analysis of publicly available HTS data, is available as example in the ``yaml_files`` folder.

Briefly, the workflow performs the following steps:

1. Quality control of Illumina reads (FastQC)
2. Trimming of the reads (e.g., removal of adapter and/or low quality sequences) (Trim galore)
3. Mapping to reference genome (BWA-MEM)
4. Convertion of mapped reads from SAM (Sequence Alignment Map) to BAM (Binary Alignment Map) format (samtools)
5. Sorting mapped reads based on read names (samtools)
6. Adding information regarding paired end reads (e.g., CIGAR field information) (samtools)
7. Re-sorting mapped reads based on chromosomal coordinates (samtools)
8. Adding basic Read-Group information regarding sample name, platform unit, platform (e.g., ILLUMINA), library and identifier (picard AddOrReplaceReadGroups)
9. Marking PCR and/or optical duplicate reads (picard MarkDuplicates)
10. Collection of summary statistics (samtools) 
11. Creation of indexes for coordinate-sorted BAM files to enable fast random access (samtools)
12. Splitting the reference genome into a predefined number of intervals for parallel processing (GATK SplitIntervals)

At this point the application of multi-sample workflow follows, during which multiple samples are concatenated into a single, unified VCF (Variant Calling Format) file, which contains the variant information for all samples:

13. Application of Base Quality Score Recalibration (BQSR) (GATK BaseRecalibrator and ApplyBQSR tools)
14. Variant calling in gVCF (genomic VCF) mode (-ERC GVCF) (GATK HaplotypeCaller)  
15. Merging of all genomic interval-split gVCF files for each sample (GATK MergeVCFs)
16. Generation of the unified VCF file (GATK CombineGVCFs and GenotypeGVCFs tools)
17. Separate annotation for SNP and INDEL variants, using the Variant Quality Score Recalibration (VQSR) method (GATK VariantRecalibrator and ApplyVQSR tools)
18. Variant filtering based on the information added during VQSR and/or custom filters (bcftools)
19. Normalization of INDELs (split multiallelic sites) (bcftools)
20. Annotation of the final dataset of filtered variants with genomic, population-related and/or clinical information (ANNOVAR)

## CWL wrappers/software tools used in this pipeline

| CWL.Subworkflow | Version | Software | CWL.wrapper | CWL.type | Docker.image |
|---|---|---|---|---|---|
| - | - | - | get-raw-files.cwl | ExpressionTool | - |
| - | - | - | split-single-paired_v2.cwl | ExpressionTool | - |
| - | 0.4.4 | Trim galore | trim-galore.cwl | CommandLineTool | kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7 |
| - | 0.11.5 | FastQC | fastqc.cwl | CommandLineTool | biowardrobe2/fastqc:v0.11.5 |
| - | latest | cp | cp.cwl | CommandLineTool | ubuntu:latest |
| - | latest | cp | rename.cwl | CommandLineTool | ubuntu:latest |
| - | - | - | check_trimming.cwl | ExpressionTool | - |
| - | latest | zcat/grep/head/cut/sed/awk | fastq_RG_extraction.cwl | CommandLineTool | ubuntu:latest |
| - | - | - | split-paired-read1-read2.cwl | ExpressionTool | - |
| - | 0.7.17-r1188 | bwa-mem | bwa-mem.cwl | CommandLineTool | quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8 |
| - | - | - | gather-bwa-sam-files.cwl | ExpressionTool | - |
| - | 1.14 | SAMtools | samtools-view.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| - | 1.14 | SAMtools | samtools-sort.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| - | 1.14 | SAMtools | samtools-fixmate.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| - | 2.26.7 | Picard tools | picard-AddOrReplaceReadGroups.cwl | CommandLineTool | quay.io/biocontainers/picard:2.26.7--hdfd78af_0 |
| - | 2.26.7 | Picard tools | picard-MarkDuplicates.cwl | CommandLineTool | quay.io/biocontainers/picard:2.26.7--hdfd78af_0 |
| - | 1.14 | SAMtools | samtools-flagstat.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| - | 4.3.0.0 | SplitIntervals (GATK4) | gatk-SplitIntervals.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| -  | 1.14 | SAMtools | samtools-index.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| gatk-bqsr-subworkflow_multiple.cwl | 4.3.0.0 | BaseRecalibrator (GATK4) | gatk-BaseRecalibrator.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | ApplyBQSR (GATK4) | gatk-ApplyBQSR_v2.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | HaplotypeCaller (GATK4) | gatk-HaplotypeCaller_multiple_v1.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | MergeVCFs (GATK4) | gatk-MergeVCFs.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | CombineGVCFs (GATK4) | gatk-CombineGVCFs.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | GenotypeGVCFs (GATK4) | gatk-GenotypeGVCFs.cwl  | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | MakeSitesOnlyVcf (GATK4) | gatk-MakeSitesOnlyVcf.cwl  | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | VariantRecalibrator (GATK4) | gatk-VariantRecalibrator.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | ApplyVQSR (GATK4) | gatk-ApplyVQSR.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 4.3.0.0 | MergeVCFs (GATK4) | gatk-MergeVCFs_vqsr.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | 1.5 | BCFtools | bcftools-view.cwl | CommandLineTool | biocontainers/bcftools:v1.5_cv3 |
| - | 1.5 | BCFtools | bcftools-norm.cwl | CommandLineTool | biocontainers/bcftools:v1.5_cv3 |
| - | 1.5 | BCFtools | bcftools-view.cwl | CommandLineTool | biocontainers/bcftools:v1.5_cv3 |
| - | 1.5 | BCFtools | bcftools-norm.cwl | CommandLineTool | biocontainers/bcftools:v1.5_cv3 |
| - | 2018-04-16 00:47:49 -0400 (Mon, 16 Apr 2018) | ANNOVAR | table-annovar.cwl | CommandLineTool | bioinfochrustrasbourg/annovar:2018Apr16 |

## Workflow structure

A tree-based representation is available below for the inspection of all workflow steps. This tree-based hierarchical structure was produced using the [data.tree](https://cran.r-project.org/web/packages/data.tree/index.html) package.

<details>
    <summary>Workflow structure</summary>

    ```bash
1   CWL-based germline variant calling (multi-sample) workflow                                                                 
2    ¦--get_raw_files                                                                                                          
3    ¦   ¦--run:                                                                                                               
4    ¦   ¦   °--../wrappers/get-raw-files.cwl                                                                                  
5    ¦   ¦--in:                                                                                                                
6    ¦   ¦   °--DIRECTORY: raw_files_directory                                                                                 
7    ¦   °--out:                                                                                                               
8    ¦       °--raw_files                                                                                                      
9    ¦--split_single_paired                                                                                                    
10   ¦   ¦--run:                                                                                                               
11   ¦   ¦   °--../wrappers/split-single-paired_v2.cwl                                                                         
12   ¦   ¦--in:                                                                                                                
13   ¦   ¦   ¦--qc_check: input_qc_check                                                                                       
14   ¦   ¦   ¦--trimming_check: input_trimming_check                                                                           
15   ¦   ¦   ¦--input_files: get_raw_files/raw_files                                                                           
16   ¦   ¦   ¦--file_split: input_file_split                                                                                   
17   ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                                             
18   ¦   ¦   °--file_split_rev: input_file_split_rev                                                                           
19   ¦   °--out:                                                                                                               
20   ¦       ¦--single_files                                                                                                   
21   ¦       ¦--paired_files                                                                                                   
22   ¦       ¦--trim_galore_for_single                                                                                         
23   ¦       ¦--trim_galore_for_paired                                                                                         
24   ¦       ¦--fastqc_for_raw                                                                                                 
25   ¦       ¦--fastqc_for_single                                                                                              
26   ¦       ¦--fastqc_for_paired                                                                                              
27   ¦       ¦--cp_command_raw                                                                                                 
28   ¦       ¦--cp_command_single                                                                                              
29   ¦       °--cp_command_paired                                                                                              
30   ¦--trim_galore_single                                                                                                     
31   ¦   ¦--run:                                                                                                               
32   ¦   ¦   °--../wrappers/trim-galore.cwl                                                                                    
33   ¦   ¦--in:                                                                                                                
34   ¦   ¦   ¦--command: split_single_paired/trim_galore_for_single                                                            
35   ¦   ¦   ¦--fq_files: split_single_paired/single_files                                                                     
36   ¦   ¦   ¦--length: tg_length                                                                                              
37   ¦   ¦   ¦--quality: tg_quality                                                                                            
38   ¦   ¦   ¦--compression: tg_compression                                                                                    
39   ¦   ¦   ¦--do_not_compress: tg_do_not_compress                                                                            
40   ¦   ¦   ¦--trim_suffix: tg_trim_suffix                                                                                    
41   ¦   ¦   ¦--strigency: tg_strigency                                                                                        
42   ¦   ¦   °--paired: $( false )                                                                                             
43   ¦   °--out:                                                                                                               
44   ¦       ¦--trim_galore                                                                                                    
45   ¦       °--trim_galore_report                                                                                             
46   ¦--trim_galore_paired                                                                                                     
47   ¦   ¦--run:                                                                                                               
48   ¦   ¦   °--../wrappers/trim-galore.cwl                                                                                    
49   ¦   ¦--in:                                                                                                                
50   ¦   ¦   ¦--command: split_single_paired/trim_galore_for_paired                                                            
51   ¦   ¦   ¦--fq_files: split_single_paired/paired_files                                                                     
52   ¦   ¦   ¦--length: tg_length                                                                                              
53   ¦   ¦   ¦--quality: tg_quality                                                                                            
54   ¦   ¦   ¦--compression: tg_compression                                                                                    
55   ¦   ¦   ¦--do_not_compress: tg_do_not_compress                                                                            
56   ¦   ¦   ¦--trim_suffix: tg_trim_suffix                                                                                    
57   ¦   ¦   ¦--strigency: tg_strigency                                                                                        
58   ¦   ¦   °--paired: $( true )                                                                                              
59   ¦   °--out:                                                                                                               
60   ¦       ¦--trim_galore                                                                                                    
61   ¦       °--trim_galore_report                                                                                             
62   ¦--fastqc_raw                                                                                                             
63   ¦   ¦--run:                                                                                                               
64   ¦   ¦   °--../wrappers/fastqc.cwl                                                                                         
65   ¦   ¦--in:                                                                                                                
66   ¦   ¦   ¦--command: split_single_paired/fastqc_for_raw                                                                    
67   ¦   ¦   °--input_files: get_raw_files/raw_files                                                                           
68   ¦   °--out:                                                                                                               
69   ¦       ¦--html_file                                                                                                      
70   ¦       °--zipped_file                                                                                                    
71   ¦--fastqc_single_trimmed                                                                                                  
72   ¦   ¦--run:                                                                                                               
73   ¦   ¦   °--../wrappers/fastqc.cwl                                                                                         
74   ¦   ¦--in:                                                                                                                
75   ¦   ¦   ¦--command: split_single_paired/fastqc_for_single                                                                 
76   ¦   ¦   °--input_files: trim_galore_single/trim_galore                                                                    
77   ¦   °--out:                                                                                                               
78   ¦       ¦--html_file                                                                                                      
79   ¦       °--zipped_file                                                                                                    
80   ¦--fastqc_paired_trimmed                                                                                                  
81   ¦   ¦--run:                                                                                                               
82   ¦   ¦   °--../wrappers/fastqc.cwl                                                                                         
83   ¦   ¦--in:                                                                                                                
84   ¦   ¦   ¦--command: split_single_paired/fastqc_for_paired                                                                 
85   ¦   ¦   °--input_files: trim_galore_paired/trim_galore                                                                    
86   ¦   °--out:                                                                                                               
87   ¦       ¦--html_file                                                                                                      
88   ¦       °--zipped_file                                                                                                    
89   ¦--cp_fastqc_raw_zip                                                                                                      
90   ¦   ¦--run:                                                                                                               
91   ¦   ¦   °--../wrappers/cp.cwl                                                                                             
92   ¦   ¦--in:                                                                                                                
93   ¦   ¦   ¦--command: split_single_paired/cp_command_raw                                                                    
94   ¦   ¦   ¦--input_files: fastqc_raw/zipped_file                                                                            
95   ¦   ¦   °--outputdir: $( "fastqc_raw_zip" )                                                                               
96   ¦   °--out:                                                                                                               
97   ¦       °--output_dir                                                                                                     
98   ¦--cp_fastqc_single_zip                                                                                                   
99   ¦   ¦--run:                                                                                                               
100  ¦   ¦   °--../wrappers/cp.cwl                                                                                             
101  ¦   ¦--in:                                                                                                                
102  ¦   ¦   ¦--command: split_single_paired/cp_command_single                                                                 
103  ¦   ¦   ¦--input_files: fastqc_single_trimmed/zipped_file                                                                 
104  ¦   ¦   °--outputdir: $( "fastqc_single_trimmed_zip" )                                                                    
105  ¦   °--out:                                                                                                               
106  ¦       °--output_dir                                                                                                     
107  ¦--cp_fastqc_paired_zip                                                                                                   
108  ¦   ¦--run:                                                                                                               
109  ¦   ¦   °--../wrappers/cp.cwl                                                                                             
110  ¦   ¦--in:                                                                                                                
111  ¦   ¦   ¦--command: split_single_paired/cp_command_paired                                                                 
112  ¦   ¦   ¦--input_files: fastqc_paired_trimmed/zipped_file                                                                 
113  ¦   ¦   °--outputdir: $( "fastqc_paired_trimmed_zip" )                                                                    
114  ¦   °--out:                                                                                                               
115  ¦       °--output_dir                                                                                                     
116  ¦--rename_fastqc_raw_html                                                                                                 
117  ¦   ¦--run:                                                                                                               
118  ¦   ¦   °--../wrappers/rename.cwl                                                                                         
119  ¦   ¦--scatter:                                                                                                           
120  ¦   ¦   °--input_file                                                                                                     
121  ¦   ¦--in:                                                                                                                
122  ¦   ¦   ¦--input_file: fastqc_raw/html_file                                                                               
123  ¦   ¦   °--run_type: $( "fastqc_raw_" )                                                                                   
124  ¦   °--out:                                                                                                               
125  ¦       °--renamed_file                                                                                                   
126  ¦--rename_fastqc_single_html                                                                                              
127  ¦   ¦--run:                                                                                                               
128  ¦   ¦   °--../wrappers/rename.cwl                                                                                         
129  ¦   ¦--scatter:                                                                                                           
130  ¦   ¦   °--input_file                                                                                                     
131  ¦   ¦--in:                                                                                                                
132  ¦   ¦   ¦--input_file: fastqc_single_trimmed/html_file                                                                    
133  ¦   ¦   °--run_type: $( "fastqc_single_trimmed_" )                                                                        
134  ¦   °--out:                                                                                                               
135  ¦       °--renamed_file                                                                                                   
136  ¦--rename_fastqc_paired_html                                                                                              
137  ¦   ¦--run:                                                                                                               
138  ¦   ¦   °--../wrappers/rename.cwl                                                                                         
139  ¦   ¦--scatter:                                                                                                           
140  ¦   ¦   °--input_file                                                                                                     
141  ¦   ¦--in:                                                                                                                
142  ¦   ¦   ¦--input_file: fastqc_paired_trimmed/html_file                                                                    
143  ¦   ¦   °--run_type: $( "fastqc_paired_trimmed_" )                                                                        
144  ¦   °--out:                                                                                                               
145  ¦       °--renamed_file                                                                                                   
146  ¦--check_trimming                                                                                                         
147  ¦   ¦--run:                                                                                                               
148  ¦   ¦   °--../wrappers/check_trimming.cwl                                                                                 
149  ¦   ¦--in:                                                                                                                
150  ¦   ¦   ¦--trimming_check: input_trimming_check                                                                           
151  ¦   ¦   ¦--file_split: input_file_split                                                                                   
152  ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                                             
153  ¦   ¦   ¦--file_split_rev: input_file_split_rev                                                                           
154  ¦   ¦   ¦--input_single: split_single_paired/single_files                                                                 
155  ¦   ¦   ¦--input_paired: split_single_paired/paired_files                                                                 
156  ¦   ¦   ¦--trimming_single: trim_galore_single/trim_galore                                                                
157  ¦   ¦   °--trimming_paired: trim_galore_paired/trim_galore                                                                
158  ¦   °--out:                                                                                                               
159  ¦       ¦--single_files                                                                                                   
160  ¦       °--paired_files                                                                                                   
161  ¦--rg_extraction_single                                                                                                   
162  ¦   ¦--run:                                                                                                               
163  ¦   ¦   °--../wrappers/fastq_RG_extraction.cwl                                                                            
164  ¦   ¦--scatter:                                                                                                           
165  ¦   ¦   °--input                                                                                                          
166  ¦   ¦--in:                                                                                                                
167  ¦   ¦   °--input: check_trimming/single_files                                                                             
168  ¦   °--out:                                                                                                               
169  ¦       °--rg_information                                                                                                 
170  ¦--bwa_mem_single                                                                                                         
171  ¦   ¦--run:                                                                                                               
172  ¦   ¦   °--../wrappers/bwa-mem.cwl                                                                                        
173  ¦   ¦--scatter:                                                                                                           
174  ¦   ¦   °--trimmed_fq_read1                                                                                               
175  ¦   ¦--scatterMethod:                                                                                                     
176  ¦   ¦   °--dotproduct                                                                                                     
177  ¦   ¦--in:                                                                                                                
178  ¦   ¦   ¦--file_split: input_file_split                                                                                   
179  ¦   ¦   ¦--sec_shorter_split_hits: bwa_mem_sec_shorter_split_hits                                                         
180  ¦   ¦   ¦--num_threads: bwa_mem_num_threads                                                                               
181  ¦   ¦   ¦--ref_genome: reference_genome                                                                                   
182  ¦   ¦   °--trimmed_fq_read1: check_trimming/single_files                                                                  
183  ¦   °--out:                                                                                                               
184  ¦       °--output                                                                                                         
185  ¦--split_paired_read1_read2                                                                                               
186  ¦   ¦--run:                                                                                                               
187  ¦   ¦   °--../wrappers/split-paired-read1-read2.cwl                                                                       
188  ¦   ¦--in:                                                                                                                
189  ¦   ¦   ¦--paired_files: check_trimming/paired_files                                                                      
190  ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                                             
191  ¦   ¦   °--file_split_rev: input_file_split_rev                                                                           
192  ¦   °--out:                                                                                                               
193  ¦       ¦--reads_1                                                                                                        
194  ¦       °--reads_2                                                                                                        
195  ¦--rg_extraction_paired                                                                                                   
196  ¦   ¦--run:                                                                                                               
197  ¦   ¦   °--../wrappers/fastq_RG_extraction.cwl                                                                            
198  ¦   ¦--scatter:                                                                                                           
199  ¦   ¦   °--input                                                                                                          
200  ¦   ¦--in:                                                                                                                
201  ¦   ¦   °--input: split_paired_read1_read2/reads_1                                                                        
202  ¦   °--out:                                                                                                               
203  ¦       °--rg_information                                                                                                 
204  ¦--bwa_mem_paired                                                                                                         
205  ¦   ¦--run:                                                                                                               
206  ¦   ¦   °--../wrappers/bwa-mem.cwl                                                                                        
207  ¦   ¦--scatter:                                                                                                           
208  ¦   ¦   ¦--trimmed_fq_read1                                                                                               
209  ¦   ¦   °--trimmed_fq_read2                                                                                               
210  ¦   ¦--scatterMethod:                                                                                                     
211  ¦   ¦   °--dotproduct                                                                                                     
212  ¦   ¦--in:                                                                                                                
213  ¦   ¦   ¦--file_split: input_file_split                                                                                   
214  ¦   ¦   ¦--sec_shorter_split_hits: bwa_mem_sec_shorter_split_hits                                                         
215  ¦   ¦   ¦--num_threads: bwa_mem_num_threads                                                                               
216  ¦   ¦   ¦--ref_genome: reference_genome                                                                                   
217  ¦   ¦   ¦--trimmed_fq_read1: split_paired_read1_read2/reads_1                                                             
218  ¦   ¦   °--trimmed_fq_read2: split_paired_read1_read2/reads_2                                                             
219  ¦   °--out:                                                                                                               
220  ¦       °--output                                                                                                         
221  ¦--gather_bwa_sam_files                                                                                                   
222  ¦   ¦--run:                                                                                                               
223  ¦   ¦   °--../wrappers/gather-bwa-sam-files.cwl                                                                           
224  ¦   ¦--in:                                                                                                                
225  ¦   ¦   ¦--single_files: bwa_mem_single/output                                                                            
226  ¦   ¦   ¦--paired_files: bwa_mem_paired/output                                                                            
227  ¦   ¦   ¦--single_files_rg_info: rg_extraction_single/rg_information                                                      
228  ¦   ¦   °--paired_files_rg_info: rg_extraction_paired/rg_information                                                      
229  ¦   °--out:                                                                                                               
230  ¦       ¦--total_sam_files                                                                                                
231  ¦       ¦--names_rgids                                                                                                    
232  ¦       ¦--names_rgpus                                                                                                    
233  ¦       °--names_rgplbs                                                                                                   
234  ¦--samtools_view_conversion                                                                                               
235  ¦   ¦--run:                                                                                                               
236  ¦   ¦   °--../wrappers/samtools-view.cwl                                                                                  
237  ¦   ¦--scatter:                                                                                                           
238  ¦   ¦   °--input                                                                                                          
239  ¦   ¦--in:                                                                                                                
240  ¦   ¦   ¦--input: gather_bwa_sam_files/total_sam_files                                                                    
241  ¦   ¦   ¦--isbam: $( true )                                                                                               
242  ¦   ¦   ¦--collapsecigar: samtools_view_collapsecigar                                                                     
243  ¦   ¦   ¦--readsingroup: samtools_view_readsingroup                                                                       
244  ¦   ¦   ¦--uncompressed: samtools_view_uncompressed                                                                       
245  ¦   ¦   ¦--readtagtostrip: samtools_view_readtagtostrip                                                                   
246  ¦   ¦   ¦--readsquality: samtools_view_readsquality                                                                       
247  ¦   ¦   ¦--readswithbits: samtools_view_readswithbits                                                                     
248  ¦   ¦   ¦--readswithoutbits: samtools_view_readswithoutbits                                                               
249  ¦   ¦   ¦--cigar: samtools_view_cigar                                                                                     
250  ¦   ¦   ¦--iscram: samtools_view_iscram                                                                                   
251  ¦   ¦   ¦--threads: samtools_view_threads                                                                                 
252  ¦   ¦   ¦--fastcompression: samtools_view_fastcompression                                                                 
253  ¦   ¦   ¦--samheader: samtools_view_samheader                                                                             
254  ¦   ¦   ¦--count: samtools_view_count                                                                                     
255  ¦   ¦   ¦--randomseed: samtools_view_randomseed                                                                           
256  ¦   ¦   ¦--region: samtools_view_region                                                                                   
257  ¦   ¦   ¦--readsinlibrary: samtools_view_readsinlibrary                                                                   
258  ¦   ¦   ¦--output_name: $( inputs.input.basename.split(".sam")[0].concat(".bam") )                                        
259  ¦   ¦   °--target_bed_file: samtools_view_target_bed_file                                                                 
260  ¦   °--out:                                                                                                               
261  ¦       °--output                                                                                                         
262  ¦--samtools_sort_by_name                                                                                                  
263  ¦   ¦--run:                                                                                                               
264  ¦   ¦   °--../wrappers/samtools-sort.cwl                                                                                  
265  ¦   ¦--scatter:                                                                                                           
266  ¦   ¦   °--input                                                                                                          
267  ¦   ¦--in:                                                                                                                
268  ¦   ¦   ¦--compression_level: samtools_sort_compression_level                                                             
269  ¦   ¦   ¦--threads: samtools_sort_threads                                                                                 
270  ¦   ¦   ¦--memory: samtools_sort_memory                                                                                   
271  ¦   ¦   ¦--input: samtools_view_conversion/output                                                                         
272  ¦   ¦   ¦--output_name: $( inputs.input.basename.split(".bam")[0].concat(".name.sorted.bam") )                            
273  ¦   ¦   °--sort_by_name: $( true )                                                                                        
274  ¦   °--out:                                                                                                               
275  ¦       °--sorted                                                                                                         
276  ¦--samtools_fixmate                                                                                                       
277  ¦   ¦--run:                                                                                                               
278  ¦   ¦   °--../wrappers/samtools-fixmate.cwl                                                                               
279  ¦   ¦--scatter:                                                                                                           
280  ¦   ¦   °--input_file                                                                                                     
281  ¦   ¦--in:                                                                                                                
282  ¦   ¦   ¦--threads: samtools_fixmate_threads                                                                              
283  ¦   ¦   ¦--output_format: samtools_fixmate_output_format                                                                  
284  ¦   ¦   ¦--input_file: samtools_sort_by_name/sorted                                                                       
285  ¦   ¦   °--output_file_name: $( inputs.input_file.basename.split(".name.sorted.bam")[0].concat(".fixed.bam") )            
286  ¦   °--out:                                                                                                               
287  ¦       °--output                                                                                                         
288  ¦--samtools_sort                                                                                                          
289  ¦   ¦--run:                                                                                                               
290  ¦   ¦   °--../wrappers/samtools-sort.cwl                                                                                  
291  ¦   ¦--scatter:                                                                                                           
292  ¦   ¦   °--input                                                                                                          
293  ¦   ¦--in:                                                                                                                
294  ¦   ¦   ¦--compression_level: samtools_sort_compression_level                                                             
295  ¦   ¦   ¦--threads: samtools_sort_threads                                                                                 
296  ¦   ¦   ¦--memory: samtools_sort_memory                                                                                   
297  ¦   ¦   ¦--input: samtools_fixmate/output                                                                                 
298  ¦   ¦   °--output_name: $( inputs.input.basename.split(".fixed.bam")[0].concat(".fixed.sorted.bam") )                     
299  ¦   °--out:                                                                                                               
300  ¦       °--sorted                                                                                                         
301  ¦--picard_addorreplacereadgroups                                                                                          
302  ¦   ¦--run:                                                                                                               
303  ¦   ¦   °--../wrappers/picard-AddOrReplaceReadGroups.cwl                                                                  
304  ¦   ¦--scatter:                                                                                                           
305  ¦   ¦   ¦--INPUT                                                                                                          
306  ¦   ¦   ¦--rgid                                                                                                           
307  ¦   ¦   ¦--rgpu                                                                                                           
308  ¦   ¦   °--rglb                                                                                                           
309  ¦   ¦--scatterMethod:                                                                                                     
310  ¦   ¦   °--dotproduct                                                                                                     
311  ¦   ¦--in:                                                                                                                
312  ¦   ¦   ¦--INPUT: samtools_sort/sorted                                                                                    
313  ¦   ¦   ¦--OUTPUT: $( inputs.INPUT.basename.split(".fixed.sorted.bam")[0].concat(".fixed.sorted.uniq.rg.bam") )           
314  ¦   ¦   ¦--rgid: gather_bwa_sam_files/names_rgids                                                                         
315  ¦   ¦   ¦--rglb: gather_bwa_sam_files/names_rgplbs                                                                        
316  ¦   ¦   ¦--rgpl: picard_addorreplacereadgroups_rgpl                                                                       
317  ¦   ¦   ¦--rgpu: gather_bwa_sam_files/names_rgpus                                                                         
318  ¦   ¦   °--rgsm: $( inputs.INPUT.basename.split(".fixed.sorted.bam")[0] )                                                 
319  ¦   °--out:                                                                                                               
320  ¦       °--output                                                                                                         
321  ¦--picard_markduplicates                                                                                                  
322  ¦   ¦--run:                                                                                                               
323  ¦   ¦   °--../wrappers/picard-MarkDuplicates.cwl                                                                          
324  ¦   ¦--scatter:                                                                                                           
325  ¦   ¦   °--INPUT                                                                                                          
326  ¦   ¦--in:                                                                                                                
327  ¦   ¦   ¦--INPUT: picard_addorreplacereadgroups/output                                                                    
328  ¦   ¦   ¦--OUTPUT: $( inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.bam")[0].concat(".fixed.sorted.uniq.rg.md.bam") )
329  ¦   ¦   °--METRICS_FILE: $( inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.bam")[0].concat(".md_metrics.txt") )       
330  ¦   °--out:                                                                                                               
331  ¦       ¦--output                                                                                                         
332  ¦       °--output_report                                                                                                  
333  ¦--samtools_flagstat                                                                                                      
334  ¦   ¦--run:                                                                                                               
335  ¦   ¦   °--../wrappers/samtools-flagstat.cwl                                                                              
336  ¦   ¦--scatter:                                                                                                           
337  ¦   ¦   °--input                                                                                                          
338  ¦   ¦--in:                                                                                                                
339  ¦   ¦   ¦--threads: samtools_flagstat_threads                                                                             
340  ¦   ¦   ¦--input: picard_markduplicates/output                                                                            
341  ¦   ¦   °--output_name: $( inputs.input.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".align.stats.txt") )    
342  ¦   °--out:                                                                                                               
343  ¦       °--output                                                                                                         
344  ¦--samtools_view_count_total                                                                                              
345  ¦   ¦--run:                                                                                                               
346  ¦   ¦   °--../wrappers/samtools-view.cwl                                                                                  
347  ¦   ¦--scatter:                                                                                                           
348  ¦   ¦   °--input                                                                                                          
349  ¦   ¦--in:                                                                                                                
350  ¦   ¦   ¦--isbam: $( false )                                                                                              
351  ¦   ¦   ¦--collapsecigar: samtools_view_collapsecigar                                                                     
352  ¦   ¦   ¦--readsingroup: samtools_view_readsingroup                                                                       
353  ¦   ¦   ¦--uncompressed: samtools_view_uncompressed                                                                       
354  ¦   ¦   ¦--readtagtostrip: samtools_view_readtagtostrip                                                                   
355  ¦   ¦   ¦--input: picard_markduplicates/output                                                                            
356  ¦   ¦   ¦--readsquality: samtools_view_readsquality                                                                       
357  ¦   ¦   ¦--cigar: samtools_view_cigar                                                                                     
358  ¦   ¦   ¦--iscram: samtools_view_iscram                                                                                   
359  ¦   ¦   ¦--threads: samtools_view_threads                                                                                 
360  ¦   ¦   ¦--fastcompression: samtools_view_fastcompression                                                                 
361  ¦   ¦   ¦--samheader: samtools_view_samheader                                                                             
362  ¦   ¦   ¦--count: TRUE                                                                                                    
363  ¦   ¦   ¦--randomseed: samtools_view_randomseed                                                                           
364  ¦   ¦   ¦--region: samtools_view_region                                                                                   
365  ¦   ¦   ¦--readsinlibrary: samtools_view_readsinlibrary                                                                   
366  ¦   ¦   °--output_name: $( inputs.input.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".count.total.txt") )    
367  ¦   °--out:                                                                                                               
368  ¦       °--output                                                                                                         
369  ¦--gatk_splitintervals                                                                                                    
370  ¦   ¦--run:                                                                                                               
371  ¦   ¦   °--../wrappers/gatk-SplitIntervals.cwl                                                                            
372  ¦   ¦--in:                                                                                                                
373  ¦   ¦   ¦--reference: reference_genome                                                                                    
374  ¦   ¦   ¦--include_intervalList: gatk_splitintervals_include_intervalList                                                 
375  ¦   ¦   ¦--exclude_intervalList: gatk_splitintervals_exclude_intervalList                                                 
376  ¦   ¦   °--scatter_count: gatk_splitintervals_scatter_count                                                               
377  ¦   °--out:                                                                                                               
378  ¦       °--intervalfiles                                                                                                  
379  ¦--samtools_index                                                                                                         
380  ¦   ¦--run:                                                                                                               
381  ¦   ¦   °--../wrappers/samtools-index.cwl                                                                                 
382  ¦   ¦--scatter:                                                                                                           
383  ¦   ¦   °--alignments                                                                                                     
384  ¦   ¦--in:                                                                                                                
385  ¦   ¦   °--alignments: picard_markduplicates/output                                                                       
386  ¦   °--out:                                                                                                               
387  ¦       °--alignments_with_index                                                                                          
388  ¦--gatk_bqsr_subworkflow                                                                                                  
389  ¦   ¦--run:                                                                                                               
390  ¦   ¦   °--../workflows/gatk-bqsr-subworkflow_multiple.cwl                                                                
391  ¦   ¦--scatter:                                                                                                           
392  ¦   ¦   °--bqsr_INPUT                                                                                                     
393  ¦   ¦--in:                                                                                                                
394  ¦   ¦   ¦--bqsr_hc_java_options: sub_bqsr_hc_java_options                                                                 
395  ¦   ¦   ¦--bqsr_INPUT: samtools_index/alignments_with_index                                                               
396  ¦   ¦   ¦--bqsr_OUTPUT: $( inputs.bqsr_INPUT.basename.split(".fixed.sorted.uniq.rg.md.bam")[0] )                          
397  ¦   ¦   ¦--bqsr_reference: reference_genome                                                                               
398  ¦   ¦   ¦--bqsr_known_sites_1: sub_bqsr_known_sites_1                                                                     
399  ¦   ¦   ¦--bqsr_known_sites_2: sub_bqsr_known_sites_2                                                                     
400  ¦   ¦   ¦--bqsr_known_sites_3: sub_bqsr_known_sites_3                                                                     
401  ¦   ¦   ¦--bqsr_intervals: gatk_splitintervals/intervalfiles                                                              
402  ¦   ¦   ¦--bqsr_interval_padding: sub_bqsr_interval_padding                                                               
403  ¦   ¦   °--bqsr_hc_native_pairHMM_threads: sub_bqsr_hc_native_pairHMM_threads                                             
404  ¦   °--out:                                                                                                               
405  ¦       ¦--o_gatk_baserecalibrator_table                                                                                  
406  ¦       ¦--o_gatk_applybqsr                                                                                               
407  ¦       ¦--o_gatk_HaplotypeCaller                                                                                         
408  ¦       °--o_gatk_MergeVCFs                                                                                               
409  ¦--gatk_CombineGVCFs                                                                                                      
410  ¦   ¦--run:                                                                                                               
411  ¦   ¦   °--../wrappers/gatk-CombineGVCFs.cwl                                                                              
412  ¦   ¦--in:                                                                                                                
413  ¦   ¦   ¦--reference: reference_genome                                                                                    
414  ¦   ¦   ¦--gvcf_files: gatk_bqsr_subworkflow/o_gatk_MergeVCFs                                                             
415  ¦   ¦   °--output_name: $( "cohort.g.vcf.gz" )                                                                            
416  ¦   °--out:                                                                                                               
417  ¦       °--output                                                                                                         
418  ¦--gatk_GenotypeGVCFs                                                                                                     
419  ¦   ¦--run:                                                                                                               
420  ¦   ¦   °--../wrappers/gatk-GenotypeGVCFs.cwl                                                                             
421  ¦   ¦--in:                                                                                                                
422  ¦   ¦   ¦--reference: reference_genome                                                                                    
423  ¦   ¦   ¦--gvcf_input: gatk_CombineGVCFs/output                                                                           
424  ¦   ¦   °--vcf_output: $( "cohort.genotyped.vcf.gz" )                                                                     
425  ¦   °--out:                                                                                                               
426  ¦       °--output                                                                                                         
427  ¦--gatk_MakeSitesOnlyVcf                                                                                                  
428  ¦   ¦--run:                                                                                                               
429  ¦   ¦   °--../wrappers/gatk-MakeSitesOnlyVcf.cwl                                                                          
430  ¦   ¦--in:                                                                                                                
431  ¦   ¦   ¦--INPUT: gatk_GenotypeGVCFs/output                                                                               
432  ¦   ¦   °--OUTPUT: $( "cohort.genotyped.sitesonly.vcf.gz" )                                                               
433  ¦   °--out:                                                                                                               
434  ¦       °--output                                                                                                         
435  ¦--gatk_VariantRecalibrator_indel                                                                                         
436  ¦   ¦--run:                                                                                                               
437  ¦   ¦   °--../wrappers/gatk-VariantRecalibrator.cwl                                                                       
438  ¦   ¦--scatter:                                                                                                           
439  ¦   ¦   °--intervals                                                                                                      
440  ¦   ¦--in:                                                                                                                
441  ¦   ¦   ¦--truth_sensitivity_trance: VariantRecalibrator_truth_sensitivity_trance_indels                                  
442  ¦   ¦   ¦--use_annotation: VariantRecalibrator_use_annotation                                                             
443  ¦   ¦   ¦--reference: reference_genome                                                                                    
444  ¦   ¦   ¦--variant: gatk_MakeSitesOnlyVcf/output                                                                          
445  ¦   ¦   ¦--mode: $( "INDEL" )                                                                                             
446  ¦   ¦   ¦--arguments_1: vqsr_arguments_indels_1                                                                           
447  ¦   ¦   ¦--resource_1: vqsr_known_sites_indels_1                                                                          
448  ¦   ¦   ¦--arguments_2: vqsr_arguments_indels_2                                                                           
449  ¦   ¦   ¦--resource_2: vqsr_known_sites_indels_2                                                                          
450  ¦   ¦   ¦--arguments_3: vqsr_arguments_indels_3                                                                           
451  ¦   ¦   ¦--resource_3: vqsr_known_sites_indels_3                                                                          
452  ¦   ¦   ¦--intervals: gatk_splitintervals/intervalfiles                                                                   
453  ¦   ¦   ¦--trust_all_polymorphic: VariantRecalibrator_trust_all_polymorphic                                               
454  ¦   ¦   °--output_name: $( "cohort.genotyped.indel" )                                                                     
455  ¦   °--out:                                                                                                               
456  ¦       ¦--recal_table                                                                                                    
457  ¦       °--tranches_file                                                                                                  
458  ¦--gatk_VariantRecalibrator_snp                                                                                           
459  ¦   ¦--run:                                                                                                               
460  ¦   ¦   °--../wrappers/gatk-VariantRecalibrator.cwl                                                                       
461  ¦   ¦--scatter:                                                                                                           
462  ¦   ¦   °--intervals                                                                                                      
463  ¦   ¦--in:                                                                                                                
464  ¦   ¦   ¦--truth_sensitivity_trance: VariantRecalibrator_truth_sensitivity_trance_snps                                    
465  ¦   ¦   ¦--use_annotation: VariantRecalibrator_use_annotation                                                             
466  ¦   ¦   ¦--reference: reference_genome                                                                                    
467  ¦   ¦   ¦--variant: gatk_MakeSitesOnlyVcf/output                                                                          
468  ¦   ¦   ¦--mode: $( "SNP" )                                                                                               
469  ¦   ¦   ¦--arguments_1: vqsr_arguments_snps_1                                                                             
470  ¦   ¦   ¦--resource_1: vqsr_known_sites_snps_1                                                                            
471  ¦   ¦   ¦--arguments_2: vqsr_arguments_snps_2                                                                             
472  ¦   ¦   ¦--resource_2: vqsr_known_sites_snps_2                                                                            
473  ¦   ¦   ¦--arguments_3: vqsr_arguments_snps_3                                                                             
474  ¦   ¦   ¦--resource_3: vqsr_known_sites_snps_3                                                                            
475  ¦   ¦   ¦--arguments_4: vqsr_arguments_snps_4                                                                             
476  ¦   ¦   ¦--resource_4: vqsr_known_sites_snps_4                                                                            
477  ¦   ¦   ¦--intervals: gatk_splitintervals/intervalfiles                                                                   
478  ¦   ¦   ¦--trust_all_polymorphic: VariantRecalibrator_trust_all_polymorphic                                               
479  ¦   ¦   °--output_name: $( "cohort.genotyped.snp" )                                                                       
480  ¦   °--out:                                                                                                               
481  ¦       ¦--recal_table                                                                                                    
482  ¦       °--tranches_file                                                                                                  
483  ¦--gatk_ApplyVQSR_indel                                                                                                   
484  ¦   ¦--run:                                                                                                               
485  ¦   ¦   °--../wrappers/gatk-ApplyVQSR.cwl                                                                                 
486  ¦   ¦--scatter:                                                                                                           
487  ¦   ¦   ¦--recal_file                                                                                                     
488  ¦   ¦   ¦--tranches_file                                                                                                  
489  ¦   ¦   °--intervals                                                                                                      
490  ¦   ¦--scatterMethod:                                                                                                     
491  ¦   ¦   °--dotproduct                                                                                                     
492  ¦   ¦--in:                                                                                                                
493  ¦   ¦   ¦--reference: reference_genome                                                                                    
494  ¦   ¦   ¦--variant: gatk_GenotypeGVCFs/output                                                                             
495  ¦   ¦   ¦--ts_filter_level: ApplyVQSR_ts_filter_level                                                                     
496  ¦   ¦   ¦--recal_file: gatk_VariantRecalibrator_indel/recal_table                                                         
497  ¦   ¦   ¦--tranches_file: gatk_VariantRecalibrator_indel/tranches_file                                                    
498  ¦   ¦   ¦--intervals: gatk_splitintervals/intervalfiles                                                                   
499  ¦   ¦   ¦--mode: $( "INDEL" )                                                                                             
500  ¦   ¦   °--output_name: $( "indel.recalibrated" )                                                                         
501  ¦   °--out:                                                                                                               
502  ¦       °--output                                                                                                         
503  ¦--gatk_ApplyVQSR_snp                                                                                                     
504  ¦   ¦--run:                                                                                                               
505  ¦   ¦   °--../wrappers/gatk-ApplyVQSR.cwl                                                                                 
506  ¦   ¦--scatter:                                                                                                           
507  ¦   ¦   ¦--variant                                                                                                        
508  ¦   ¦   ¦--recal_file                                                                                                     
509  ¦   ¦   ¦--tranches_file                                                                                                  
510  ¦   ¦   °--intervals                                                                                                      
511  ¦   ¦--scatterMethod:                                                                                                     
512  ¦   ¦   °--dotproduct                                                                                                     
513  ¦   ¦--in:                                                                                                                
514  ¦   ¦   ¦--reference: reference_genome                                                                                    
515  ¦   ¦   ¦--variant: gatk_ApplyVQSR_indel/output                                                                           
516  ¦   ¦   ¦--ts_filter_level: ApplyVQSR_ts_filter_level                                                                     
517  ¦   ¦   ¦--recal_file: gatk_VariantRecalibrator_snp/recal_table                                                           
518  ¦   ¦   ¦--tranches_file: gatk_VariantRecalibrator_snp/tranches_file                                                      
519  ¦   ¦   ¦--intervals: gatk_splitintervals/intervalfiles                                                                   
520  ¦   ¦   ¦--mode: $( "SNP" )                                                                                               
521  ¦   ¦   °--output_name: $( "snp.recalibrated" )                                                                           
522  ¦   °--out:                                                                                                               
523  ¦       °--output                                                                                                         
524  ¦--gatk_VQSR_MergeVCFs                                                                                                    
525  ¦   ¦--run:                                                                                                               
526  ¦   ¦   °--../wrappers/gatk-MergeVCFs_vqsr.cwl                                                                            
527  ¦   ¦--in:                                                                                                                
528  ¦   ¦   ¦--INPUT: gatk_ApplyVQSR_snp/output                                                                               
529  ¦   ¦   °--OUTPUT: $( "cohort.snp.indel.vqsr.recalibrated" )                                                              
530  ¦   °--out:                                                                                                               
531  ¦       °--output                                                                                                         
532  ¦--bcftools_view_filter_vqsr                                                                                              
533  ¦   ¦--run:                                                                                                               
534  ¦   ¦   °--../wrappers/bcftools-view.cwl                                                                                  
535  ¦   ¦--in:                                                                                                                
536  ¦   ¦   ¦--input: gatk_VQSR_MergeVCFs/output                                                                              
537  ¦   ¦   ¦--threads: bcftools_view_threads                                                                                 
538  ¦   ¦   ¦--include: bcftools_view_include_VQSR_filters                                                                    
539  ¦   ¦   ¦--output_type: $( "z" )                                                                                          
540  ¦   ¦   °--output_name: $("cohort.snp.indel.vqsr.recalibrated.filtered.vcf.gz")                                           
541  ¦   °--out:                                                                                                               
542  ¦       °--output                                                                                                         
543  ¦--bcftools_norm_vqsr                                                                                                     
544  ¦   ¦--run:                                                                                                               
545  ¦   ¦   °--../wrappers/bcftools-norm.cwl                                                                                  
546  ¦   ¦--in:                                                                                                                
547  ¦   ¦   ¦--input: gatk_VQSR_MergeVCFs/output                                                                              
548  ¦   ¦   ¦--threads: bcftools_norm_threads                                                                                 
549  ¦   ¦   ¦--reference: reference_genome                                                                                    
550  ¦   ¦   ¦--multiallelics: bcftools_norm_multiallelics                                                                     
551  ¦   ¦   °--output_type: $( "v" )                                                                                          
552  ¦   °--out:                                                                                                               
553  ¦       °--output                                                                                                         
554  °--table_annovar_filtered                                                                                                 
555      ¦--run:                                                                                                               
556      ¦   °--../wrappers/table-annovar.cwl                                                                                  
557      ¦--in:                                                                                                                
558      ¦   ¦--query_file: bcftools_norm_vqsr/output                                                                          
559      ¦   ¦--database_location: table_annovar_database_location                                                             
560      ¦   ¦--build_over: table_annovar_build_over                                                                           
561      ¦   ¦--remove: table_annovar_remove                                                                                   
562      ¦   ¦--protocol: table_annovar_protocol                                                                               
563      ¦   ¦--operation: table_annovar_operation                                                                             
564      ¦   ¦--na_string: table_annovar_na_string                                                                             
565      ¦   ¦--vcfinput: table_annovar_vcfinput                                                                               
566      ¦   ¦--otherinfo: table_annovar_otherinfo                                                                             
567      ¦   °--convert_arg: table_annovar_convert_arg                                                                         
568      °--out:                                                                                                               
569          ¦--multianno_vcf                                                                                                  
570          ¦--multianno_txt                                                                                                  
571          °--avinput                                                                                                        
```

</details>