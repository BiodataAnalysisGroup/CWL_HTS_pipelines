# CWL-based (single-sample) workflow for germline variant calling

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

At this point the application of single-sample workflow follows, during which multiple samples are accepted as input and they are not merged into a unified VCF file but are rather processed separately in each step of the workflow, leading to the production of a VCF file for each sample:

13. Application of Base Quality Score Recalibration (BQSR) (GATK BaseRecalibrator, GatherBQSRReports and ApplyBQSR tools)
14. Variant calling (GATK HaplotypeCaller)  
15. Merging of all genomic interval-split gVCF files for each sample (GATK MergeVCFs)
16. Separate annotation of SNPs and INDELs based on pretrained Convolutional Neural Network (CNN) models (GATK SelectVariants, CNNScoreVariants and FilterVariantTranches tools)
17. (Optional) Independent step of hard-filtering (GATK VariantFiltration)
18. Variant filtering based on the information added during VQSR and/or custom filters (bcftools)
19. Normalization of INDELs (split multiallelic sites) (bcftools)
20. Annotation of the final dataset of filtered variants with genomic, population-related and/or clinical information (ANNOVAR)

## CWL wrappers/software tools used in this pipeline

| CWL.Subworkflow | Software | CWL.wrapper | CWL.type | Docker.image |
|---|---|---|---|---|
| - | - | get-raw-files.cwl | ExpressionTool | - |
| - | - | split-single-paired_v2.cwl | ExpressionTool | - |
| - | Trim galore | trim-galore.cwl | CommandLineTool | kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7 |
| - | FastQC | fastqc.cwl | CommandLineTool | biowardrobe2/fastqc:v0.11.5 |
| - | cp | cp.cwl | CommandLineTool | ubuntu:latest |
| - | cp | rename.cwl | CommandLineTool | ubuntu:latest |
| - | - | check_trimming.cwl | ExpressionTool | - |
| - | zcat/grep/head/cut/sed/awk | fastq_RG_extraction.cwl | CommandLineTool | ubuntu:latest |
| - | - | split-paired-read1-read2.cwl | ExpressionTool | - |
| - | bwa-mem | bwa-mem.cwl | CommandLineTool | quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8 |
| - | - | gather-bwa-sam-files_single_v2.cwl | ExpressionTool | - |
| - | SAMtools | samtools-view.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| - | SAMtools | samtools-sort.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| - | SAMtools | samtools-fixmate.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| - | Picard tools | picard-AddOrReplaceReadGroups.cwl | CommandLineTool | quay.io/biocontainers/picard:2.26.7--hdfd78af_0 |
| - | Picard tools | picard-MarkDuplicates.cwl | CommandLineTool | quay.io/biocontainers/picard:2.26.7--hdfd78af_0 |
| - | SAMtools | samtools-flagstat.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| - | GATK4 | gatk-SplitIntervals.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | SAMtools | samtools-index.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| gatk-bqsr-subworkflow.cwl | GATK4 | gatk-BaseRecalibrator.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | GATK4 | gatk-GatherBQSRReports.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | GATK4 | gatk-ApplyBQSR.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| gatk-haplotypecaller-subworkflow.cwl | GATK4 | gatk-HaplotypeCaller_single_v1.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | GATK4 | gatk-MergeVCFs.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | GATK4 | gatk-SelectVariants.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | GATK4 | gatk-VariantFiltration.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | bgzip | bgzip.cwl | CommandLineTool | bioslimcontainers/tabix:1.7 |
| - | tabix | tabix.cwl | CommandLineTool | bioslimcontainers/tabix:1.7 |
| - | BCFtools | bcftools-concat.cwl | CommandLineTool | biocontainers/bcftools:v1.5_cv3 |
| - | BCFtools | bcftools-view.cwl | CommandLineTool | biocontainers/bcftools:v1.5_cv3 |
| - | BCFtools | bcftools-norm.cwl | CommandLineTool | biocontainers/bcftools:v1.5_cv3 |
| - | GATK4 | gatk-CNNScoreVariants.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | GATK4 | gatk-FilterVariantTranches.cwl | CommandLineTool | broadinstitute/gatk:4.3.0.0 |
| - | ANNOVAR | table-annovar.cwl | CommandLineTool | bioinfochrustrasbourg/annovar:2018Apr16 |

## Workflow structure

A tree-based representation is available below for the inspection of all workflow steps. This tree-based hierarchical structure was produced using the [data.tree](https://cran.r-project.org/web/packages/data.tree/index.html) package.

<details>
    <summary>Workflow structure</summary>

    ```bash
1   CWL-based germline variant calling (single-sample) workflow                                                              
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
175  ¦   ¦--in:                                                                                                              
176  ¦   ¦   ¦--file_split: input_file_split                                                                                 
177  ¦   ¦   ¦--sec_shorter_split_hits: bwa_mem_sec_shorter_split_hits                                                       
178  ¦   ¦   ¦--num_threads: bwa_mem_num_threads                                                                             
179  ¦   ¦   ¦--ref_genome: reference_genome                                                                                 
180  ¦   ¦   °--trimmed_fq_read1: check_trimming/single_files                                                                
181  ¦   °--out:                                                                                                             
182  ¦       °--output                                                                                                       
183  ¦--split_paired_read1_read2                                                                                             
184  ¦   ¦--run:                                                                                                             
185  ¦   ¦   °--../wrappers/split-paired-read1-read2.cwl                                                                     
186  ¦   ¦--in:                                                                                                              
187  ¦   ¦   ¦--paired_files: check_trimming/paired_files                                                                    
188  ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                                           
189  ¦   ¦   °--file_split_rev: input_file_split_rev                                                                         
190  ¦   °--out:                                                                                                             
191  ¦       ¦--reads_1                                                                                                      
192  ¦       °--reads_2                                                                                                      
193  ¦--rg_extraction_paired                                                                                                 
194  ¦   ¦--run:                                                                                                             
195  ¦   ¦   °--../wrappers/fastq_RG_extraction.cwl                                                                          
196  ¦   ¦--scatter:                                                                                                         
197  ¦   ¦   °--input                                                                                                        
198  ¦   ¦--in:                                                                                                              
199  ¦   ¦   °--input: split_paired_read1_read2/reads_1                                                                      
200  ¦   °--out:                                                                                                             
201  ¦       °--rg_information                                                                                               
202  ¦--bwa_mem_paired                                                                                                       
203  ¦   ¦--run:                                                                                                             
204  ¦   ¦   °--../wrappers/bwa-mem.cwl                                                                                      
205  ¦   ¦--scatter:                                                                                                         
206  ¦   ¦   ¦--trimmed_fq_read1                                                                                             
207  ¦   ¦   °--trimmed_fq_read2                                                                                             
208  ¦   ¦--scatterMethod:                                                                                                   
209  ¦   ¦   °--dotproduct                                                                                                   
210  ¦   ¦--in:                                                                                                              
211  ¦   ¦   ¦--file_split: input_file_split                                                                                 
212  ¦   ¦   ¦--sec_shorter_split_hits: bwa_mem_sec_shorter_split_hits                                                       
213  ¦   ¦   ¦--num_threads: bwa_mem_num_threads                                                                             
214  ¦   ¦   ¦--ref_genome: reference_genome                                                                                 
215  ¦   ¦   ¦--trimmed_fq_read1: split_paired_read1_read2/reads_1                                                           
216  ¦   ¦   °--trimmed_fq_read2: split_paired_read1_read2/reads_2                                                           
217  ¦   °--out:                                                                                                             
218  ¦       °--output                                                                                                       
219  ¦--gather_bwa_sam_files                                                                                                 
220  ¦   ¦--run:                                                                                                             
221  ¦   ¦   °--../wrappers/gather-bwa-sam-files_single_v2.cwl                                                               
222  ¦   ¦--in:                                                                                                              
223  ¦   ¦   ¦--single_files: bwa_mem_single/output                                                                          
224  ¦   ¦   ¦--paired_files: bwa_mem_paired/output                                                                          
225  ¦   ¦   ¦--single_files_rg_info: rg_extraction_single/rg_information                                                    
226  ¦   ¦   °--paired_files_rg_info: rg_extraction_paired/rg_information                                                    
227  ¦   °--out:                                                                                                             
228  ¦       ¦--total_sam_files                                                                                              
229  ¦       ¦--names_rgids                                                                                                  
230  ¦       ¦--names_rgpus                                                                                                  
231  ¦       °--names_rgplbs                                                                                                 
232  ¦--samtools_view_conversion                                                                                             
233  ¦   ¦--run:                                                                                                             
234  ¦   ¦   °--../wrappers/samtools-view.cwl                                                                                
235  ¦   ¦--scatter:                                                                                                         
236  ¦   ¦   °--input                                                                                                        
237  ¦   ¦--in:                                                                                                              
238  ¦   ¦   ¦--isbam: $( true )                                                                                             
239  ¦   ¦   ¦--collapsecigar: samtools_view_collapsecigar                                                                   
240  ¦   ¦   ¦--readsingroup: samtools_view_readsingroup                                                                     
241  ¦   ¦   ¦--uncompressed: samtools_view_uncompressed                                                                     
242  ¦   ¦   ¦--readtagtostrip: samtools_view_readtagtostrip                                                                 
243  ¦   ¦   ¦--input: gather_bwa_sam_files/total_sam_files                                                                  
244  ¦   ¦   ¦--readsquality: samtools_view_readsquality                                                                     
245  ¦   ¦   ¦--readswithbits: samtools_view_readswithbits                                                                   
246  ¦   ¦   ¦--readswithoutbits: samtools_view_readswithoutbits                                                             
247  ¦   ¦   ¦--cigar: samtools_view_cigar                                                                                   
248  ¦   ¦   ¦--iscram: samtools_view_iscram                                                                                 
249  ¦   ¦   ¦--threads: samtools_view_threads                                                                               
250  ¦   ¦   ¦--fastcompression: samtools_view_fastcompression                                                               
251  ¦   ¦   ¦--samheader: samtools_view_samheader                                                                           
252  ¦   ¦   ¦--count: samtools_view_count                                                                                   
253  ¦   ¦   ¦--randomseed: samtools_view_randomseed                                                                         
254  ¦   ¦   ¦--region: samtools_view_region                                                                                 
255  ¦   ¦   ¦--readsinlibrary: samtools_view_readsinlibrary                                                                 
256  ¦   ¦   °--output_name: $(inputs.input.basename.split(".sam")[0].concat(".bam"))                                        
257  ¦   °--out:                                                                                                             
258  ¦       °--output                                                                                                       
259  ¦--samtools_sort_by_name                                                                                                
260  ¦   ¦--run:                                                                                                             
261  ¦   ¦   °--../wrappers/samtools-sort.cwl                                                                                
262  ¦   ¦--scatter:                                                                                                         
263  ¦   ¦   °--input                                                                                                        
264  ¦   ¦--in:                                                                                                              
265  ¦   ¦   ¦--compression_level: samtools_sort_compression_level                                                           
266  ¦   ¦   ¦--threads: samtools_sort_threads                                                                               
267  ¦   ¦   ¦--memory: samtools_sort_memory                                                                                 
268  ¦   ¦   ¦--input: samtools_view_conversion/output                                                                       
269  ¦   ¦   ¦--output_name: $( inputs.input.basename.split(".bam")[0].concat(".name.sorted.bam") )                          
270  ¦   ¦   °--sort_by_name: $( true )                                                                                      
271  ¦   °--out:                                                                                                             
272  ¦       °--sorted                                                                                                       
273  ¦--samtools_fixmate                                                                                                     
274  ¦   ¦--run:                                                                                                             
275  ¦   ¦   °--../wrappers/samtools-fixmate.cwl                                                                             
276  ¦   ¦--scatter:                                                                                                         
277  ¦   ¦   °--input_file                                                                                                   
278  ¦   ¦--in:                                                                                                              
279  ¦   ¦   ¦--threads: samtools_fixmate_threads                                                                            
280  ¦   ¦   ¦--output_format: samtools_fixmate_output_format                                                                
281  ¦   ¦   ¦--input_file: samtools_sort_by_name/sorted                                                                     
282  ¦   ¦   °--output_file_name: $(inputs.input_file.basename.split(".name.sorted.bam")[0].concat(".fixed.bam"))            
283  ¦   °--out:                                                                                                             
284  ¦       °--output                                                                                                       
285  ¦--samtools_sort                                                                                                        
286  ¦   ¦--run:                                                                                                             
287  ¦   ¦   °--../wrappers/samtools-sort.cwl                                                                                
288  ¦   ¦--scatter:                                                                                                         
289  ¦   ¦   °--input                                                                                                        
290  ¦   ¦--in:                                                                                                              
291  ¦   ¦   ¦--compression_level: samtools_sort_compression_level                                                           
292  ¦   ¦   ¦--threads: samtools_sort_threads                                                                               
293  ¦   ¦   ¦--memory: samtools_sort_memory                                                                                 
294  ¦   ¦   ¦--input: samtools_fixmate/output                                                                               
295  ¦   ¦   °--output_name: $(inputs.input.basename.split(".fixed.bam")[0].concat(".fixed.sorted.bam"))                     
296  ¦   °--out:                                                                                                             
297  ¦       °--sorted                                                                                                       
298  ¦--picard_addorreplacereadgroups                                                                                        
299  ¦   ¦--run:                                                                                                             
300  ¦   ¦   °--../wrappers/picard-AddOrReplaceReadGroups.cwl                                                                
301  ¦   ¦--scatter:                                                                                                         
302  ¦   ¦   ¦--INPUT                                                                                                        
303  ¦   ¦   ¦--rgid                                                                                                         
304  ¦   ¦   ¦--rgpu                                                                                                         
305  ¦   ¦   °--rglb                                                                                                         
306  ¦   ¦--scatterMethod:                                                                                                   
307  ¦   ¦   °--dotproduct                                                                                                   
308  ¦   ¦--in:                                                                                                              
309  ¦   ¦   ¦--INPUT: samtools_sort/sorted                                                                                  
310  ¦   ¦   ¦--OUTPUT: $(inputs.INPUT.basename.split(".fixed.sorted.bam")[0].concat(".fixed.sorted.uniq.rg.bam"))           
311  ¦   ¦   ¦--rgid: gather_bwa_sam_files/names_rgids                                                                       
312  ¦   ¦   ¦--rglb: gather_bwa_sam_files/names_rgplbs                                                                      
313  ¦   ¦   ¦--rgpl: picard_addorreplacereadgroups_rgpl                                                                     
314  ¦   ¦   ¦--rgpu: gather_bwa_sam_files/names_rgpus                                                                       
315  ¦   ¦   °--rgsm: $(inputs.INPUT.basename.split(".fixed.sorted.bam")[0])                                                 
316  ¦   °--out:                                                                                                             
317  ¦       °--output                                                                                                       
318  ¦--picard_markduplicates                                                                                                
319  ¦   ¦--run:                                                                                                             
320  ¦   ¦   °--../wrappers/picard-MarkDuplicates.cwl                                                                        
321  ¦   ¦--scatter:                                                                                                         
322  ¦   ¦   °--INPUT                                                                                                        
323  ¦   ¦--in:                                                                                                              
324  ¦   ¦   ¦--INPUT: picard_addorreplacereadgroups/output                                                                  
325  ¦   ¦   ¦--OUTPUT: $(inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.bam")[0].concat(".fixed.sorted.uniq.rg.md.bam"))
326  ¦   ¦   °--METRICS_FILE: $(inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.bam")[0].concat(".md_metrics.txt"))       
327  ¦   °--out:                                                                                                             
328  ¦       ¦--output                                                                                                       
329  ¦       °--output_report                                                                                                
330  ¦--samtools_flagstat                                                                                                    
331  ¦   ¦--run:                                                                                                             
332  ¦   ¦   °--../wrappers/samtools-flagstat.cwl                                                                            
333  ¦   ¦--scatter:                                                                                                         
334  ¦   ¦   °--input                                                                                                        
335  ¦   ¦--in:                                                                                                              
336  ¦   ¦   ¦--threads: samtools_flagstat_threads                                                                           
337  ¦   ¦   ¦--input: picard_markduplicates/output                                                                          
338  ¦   ¦   °--output_name: $(inputs.input.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".align.stats.txt"))    
339  ¦   °--out:                                                                                                             
340  ¦       °--output                                                                                                       
341  ¦--samtools_view_count_total                                                                                            
342  ¦   ¦--run:                                                                                                             
343  ¦   ¦   °--../wrappers/samtools-view.cwl                                                                                
344  ¦   ¦--scatter:                                                                                                         
345  ¦   ¦   °--input                                                                                                        
346  ¦   ¦--in:                                                                                                              
347  ¦   ¦   ¦--isbam: $( false )                                                                                            
348  ¦   ¦   ¦--collapsecigar: samtools_view_collapsecigar                                                                   
349  ¦   ¦   ¦--readsingroup: samtools_view_readsingroup                                                                     
350  ¦   ¦   ¦--uncompressed: samtools_view_uncompressed                                                                     
351  ¦   ¦   ¦--readtagtostrip: samtools_view_readtagtostrip                                                                 
352  ¦   ¦   ¦--input: picard_markduplicates/output                                                                          
353  ¦   ¦   ¦--readsquality: samtools_view_readsquality                                                                     
354  ¦   ¦   ¦--readswithbits: samtools_view_readswithbits                                                                   
355  ¦   ¦   ¦--cigar: samtools_view_cigar                                                                                   
356  ¦   ¦   ¦--iscram: samtools_view_iscram                                                                                 
357  ¦   ¦   ¦--threads: samtools_view_threads                                                                               
358  ¦   ¦   ¦--fastcompression: samtools_view_fastcompression                                                               
359  ¦   ¦   ¦--samheader: samtools_view_samheader                                                                           
360  ¦   ¦   ¦--count: TRUE                                                                                                  
361  ¦   ¦   ¦--randomseed: samtools_view_randomseed                                                                         
362  ¦   ¦   ¦--region: samtools_view_region                                                                                 
363  ¦   ¦   ¦--readsinlibrary: samtools_view_readsinlibrary                                                                 
364  ¦   ¦   °--output_name: $(inputs.input.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".count.total.txt"))    
365  ¦   °--out:                                                                                                             
366  ¦       °--output                                                                                                       
367  ¦--gatk_splitintervals                                                                                                  
368  ¦   ¦--run:                                                                                                             
369  ¦   ¦   °--../wrappers/gatk-SplitIntervals.cwl                                                                          
370  ¦   ¦--in:                                                                                                              
371  ¦   ¦   ¦--reference: reference_genome                                                                                  
372  ¦   ¦   ¦--include_intervalList: gatk_splitintervals_include_intervalList                                               
373  ¦   ¦   ¦--exclude_intervalList: gatk_splitintervals_exclude_intervalList                                               
374  ¦   ¦   °--scatter_count: gatk_splitintervals_scatter_count                                                             
375  ¦   °--out:                                                                                                             
376  ¦       °--intervalfiles                                                                                                
377  ¦--samtools_index                                                                                                       
378  ¦   ¦--run:                                                                                                             
379  ¦   ¦   °--../wrappers/samtools-index.cwl                                                                               
380  ¦   ¦--scatter:                                                                                                         
381  ¦   ¦   °--alignments                                                                                                   
382  ¦   ¦--in:                                                                                                              
383  ¦   ¦   °--alignments: picard_markduplicates/output                                                                     
384  ¦   °--out:                                                                                                             
385  ¦       °--alignments_with_index                                                                                        
386  ¦--gatk_bqsr_subworkflow                                                                                                
387  ¦   ¦--run:                                                                                                             
388  ¦   ¦   °--../workflows/gatk-bqsr-subworkflow.cwl                                                                       
389  ¦   ¦--scatter:                                                                                                         
390  ¦   ¦   °--bqsr_INPUT                                                                                                   
391  ¦   ¦--in:                                                                                                              
392  ¦   ¦   ¦--bqsr_INPUT: samtools_index/alignments_with_index                                                             
393  ¦   ¦   ¦--bqsr_OUTPUT: $(inputs.bqsr_INPUT.basename.split(".fixed.sorted.uniq.rg.md.bam")[0])                          
394  ¦   ¦   ¦--bqsr_reference: reference_genome                                                                             
395  ¦   ¦   ¦--bqsr_known_sites_1: sub_bqsr_known_sites_1                                                                   
396  ¦   ¦   ¦--bqsr_known_sites_2: sub_bqsr_known_sites_2                                                                   
397  ¦   ¦   ¦--bqsr_known_sites_3: sub_bqsr_known_sites_3                                                                   
398  ¦   ¦   ¦--bqsr_intervals: gatk_splitintervals/intervalfiles                                                            
399  ¦   ¦   °--bqsr_interval_padding: sub_bqsr_interval_padding                                                             
400  ¦   °--out:                                                                                                             
401  ¦       °--o_gatk_gatherbqsrreports                                                                                     
402  ¦--gatk_applybqsr                                                                                                       
403  ¦   ¦--run:                                                                                                             
404  ¦   ¦   °--../wrappers/gatk-ApplyBQSR.cwl                                                                               
405  ¦   ¦--scatter:                                                                                                         
406  ¦   ¦   ¦--INPUT                                                                                                        
407  ¦   ¦   °--bqsr_recal_file                                                                                              
408  ¦   ¦--scatterMethod:                                                                                                   
409  ¦   ¦   °--dotproduct                                                                                                   
410  ¦   ¦--in:                                                                                                              
411  ¦   ¦   ¦--INPUT: samtools_index/alignments_with_index                                                                  
412  ¦   ¦   ¦--bqsr_recal_file: gatk_bqsr_subworkflow/o_gatk_gatherbqsrreports                                              
413  ¦   ¦   ¦--OUTPUT: $(inputs.INPUT.basename.split(".fixed.sorted.uniq.rg.md.bam")[0].concat(".bqsr.bam"))                
414  ¦   ¦   °--reference: reference_genome                                                                                  
415  ¦   °--out:                                                                                                             
416  ¦       °--output                                                                                                       
417  ¦--samtools_index_2                                                                                                     
418  ¦   ¦--run:                                                                                                             
419  ¦   ¦   °--../wrappers/samtools-index.cwl                                                                               
420  ¦   ¦--scatter:                                                                                                         
421  ¦   ¦   °--alignments                                                                                                   
422  ¦   ¦--in:                                                                                                              
423  ¦   ¦   °--alignments: gatk_applybqsr/output                                                                            
424  ¦   °--out:                                                                                                             
425  ¦       °--alignments_with_index                                                                                        
426  ¦--gatk_haplotypecaller_subworkflow                                                                                     
427  ¦   ¦--run:                                                                                                             
428  ¦   ¦   °--../workflows/gatk-haplotypecaller-subworkflow.cwl                                                            
429  ¦   ¦--scatter:                                                                                                         
430  ¦   ¦   °--hf_INPUT                                                                                                     
431  ¦   ¦--in:                                                                                                              
432  ¦   ¦   ¦--hf_INPUT: samtools_index/alignments_with_index                                                               
433  ¦   ¦   ¦--hf_OUTPUT: $(inputs.hf_INPUT.basename.split(".fixed.sorted.uniq.rg.md.bam")[0])                              
434  ¦   ¦   ¦--hf_reference: reference_genome                                                                               
435  ¦   ¦   ¦--hf_intervals: gatk_splitintervals/intervalfiles                                                              
436  ¦   ¦   ¦--hf_hc_native_pairHMM_threads: sub_hc_native_pairHMM_threads                                                  
437  ¦   ¦   °--hf_hc_java_options: sub_hc_java_options                                                                      
438  ¦   °--out:                                                                                                             
439  ¦       °--o_gatk_MergeVCFs                                                                                             
440  ¦--gatk_SelectVariants_snps                                                                                             
441  ¦   ¦--run:                                                                                                             
442  ¦   ¦   °--../wrappers/gatk-SelectVariants.cwl                                                                          
443  ¦   ¦--scatter:                                                                                                         
444  ¦   ¦   °--variant                                                                                                      
445  ¦   ¦--in:                                                                                                              
446  ¦   ¦   ¦--reference: reference_genome                                                                                  
447  ¦   ¦   ¦--variant: gatk_haplotypecaller_subworkflow/o_gatk_MergeVCFs                                                   
448  ¦   ¦   ¦--OUTPUT: $(inputs.variant.basename.split(".vcf.gz")[0].concat(".snp.vcf"))                                    
449  ¦   ¦   °--variant_type: $( "SNP" )                                                                                     
450  ¦   °--out:                                                                                                             
451  ¦       °--output                                                                                                       
452  ¦--gatk_SelectVariants_indels                                                                                           
453  ¦   ¦--run:                                                                                                             
454  ¦   ¦   °--../wrappers/gatk-SelectVariants.cwl                                                                          
455  ¦   ¦--scatter:                                                                                                         
456  ¦   ¦   °--variant                                                                                                      
457  ¦   ¦--in:                                                                                                              
458  ¦   ¦   ¦--reference: reference_genome                                                                                  
459  ¦   ¦   ¦--variant: gatk_haplotypecaller_subworkflow/o_gatk_MergeVCFs                                                   
460  ¦   ¦   ¦--OUTPUT: $(inputs.variant.basename.split(".vcf.gz")[0].concat(".indel.vcf"))                                  
461  ¦   ¦   °--variant_type: $( "INDEL" )                                                                                   
462  ¦   °--out:                                                                                                             
463  ¦       °--output                                                                                                       
464  ¦--gatk_VariantFiltration_snps                                                                                          
465  ¦   ¦--run:                                                                                                             
466  ¦   ¦   °--../wrappers/gatk-VariantFiltration.cwl                                                                       
467  ¦   ¦--scatter:                                                                                                         
468  ¦   ¦   °--variant                                                                                                      
469  ¦   ¦--in:                                                                                                              
470  ¦   ¦   ¦--reference: reference_genome                                                                                  
471  ¦   ¦   ¦--variant: gatk_SelectVariants_snps/output                                                                     
472  ¦   ¦   ¦--window: VariantFiltration_window                                                                             
473  ¦   ¦   ¦--cluster: VariantFiltration_cluster                                                                           
474  ¦   ¦   ¦--filter_name: VariantFiltration_filter_name_snp                                                               
475  ¦   ¦   ¦--filter: VariantFiltration_filter_snp                                                                         
476  ¦   ¦   °--OUTPUT: $(inputs.variant.basename.split(".snp.vcf")[0].concat(".filtered.snp.vcf"))                          
477  ¦   °--out:                                                                                                             
478  ¦       °--output                                                                                                       
479  ¦--gatk_VariantFiltration_indels                                                                                        
480  ¦   ¦--run:                                                                                                             
481  ¦   ¦   °--../wrappers/gatk-VariantFiltration.cwl                                                                       
482  ¦   ¦--scatter:                                                                                                         
483  ¦   ¦   °--variant                                                                                                      
484  ¦   ¦--in:                                                                                                              
485  ¦   ¦   ¦--reference: reference_genome                                                                                  
486  ¦   ¦   ¦--variant: gatk_SelectVariants_indels/output                                                                   
487  ¦   ¦   ¦--window: VariantFiltration_window                                                                             
488  ¦   ¦   ¦--cluster: VariantFiltration_cluster                                                                           
489  ¦   ¦   ¦--filter_name: VariantFiltration_filter_name_indel                                                             
490  ¦   ¦   ¦--filter: VariantFiltration_filter_indel                                                                       
491  ¦   ¦   °--OUTPUT: $(inputs.variant.basename.split(".indel.vcf")[0].concat(".filtered.indel.vcf"))                      
492  ¦   °--out:                                                                                                             
493  ¦       °--output                                                                                                       
494  ¦--bgzip_snps                                                                                                           
495  ¦   ¦--run:                                                                                                             
496  ¦   ¦   °--../wrappers/bgzip.cwl                                                                                        
497  ¦   ¦--scatter:                                                                                                         
498  ¦   ¦   °--input                                                                                                        
499  ¦   ¦--in:                                                                                                              
500  ¦   ¦   °--input: gatk_VariantFiltration_snps/output                                                                    
501  ¦   °--out:                                                                                                             
502  ¦       °--output                                                                                                       
503  ¦--tabix_snps                                                                                                           
504  ¦   ¦--run:                                                                                                             
505  ¦   ¦   °--../wrappers/tabix.cwl                                                                                        
506  ¦   ¦--scatter:                                                                                                         
507  ¦   ¦   °--input                                                                                                        
508  ¦   ¦--in:                                                                                                              
509  ¦   ¦   °--input: bgzip_snps/output                                                                                     
510  ¦   °--out:                                                                                                             
511  ¦       °--output                                                                                                       
512  ¦--bgzip_indels                                                                                                         
513  ¦   ¦--run:                                                                                                             
514  ¦   ¦   °--../wrappers/bgzip.cwl                                                                                        
515  ¦   ¦--scatter:                                                                                                         
516  ¦   ¦   °--input                                                                                                        
517  ¦   ¦--in:                                                                                                              
518  ¦   ¦   °--input: gatk_VariantFiltration_snps/output                                                                    
519  ¦   °--out:                                                                                                             
520  ¦       °--output                                                                                                       
521  ¦--tabix_indels                                                                                                         
522  ¦   ¦--run:                                                                                                             
523  ¦   ¦   °--../wrappers/tabix.cwl                                                                                        
524  ¦   ¦--scatter:                                                                                                         
525  ¦   ¦   °--input                                                                                                        
526  ¦   ¦--in:                                                                                                              
527  ¦   ¦   °--input: bgzip_indels/output                                                                                   
528  ¦   °--out:                                                                                                             
529  ¦       °--output                                                                                                       
530  ¦--bcftools_concat                                                                                                      
531  ¦   ¦--run:                                                                                                             
532  ¦   ¦   °--../wrappers/bcftools-concat.cwl                                                                              
533  ¦   ¦--scatter:                                                                                                         
534  ¦   ¦   ¦--input1                                                                                                       
535  ¦   ¦   °--input2                                                                                                       
536  ¦   ¦--scatterMethod:                                                                                                   
537  ¦   ¦   °--dotproduct                                                                                                   
538  ¦   ¦--in:                                                                                                              
539  ¦   ¦   ¦--input1: tabix_snps/output                                                                                    
540  ¦   ¦   ¦--input2: tabix_indels/output                                                                                  
541  ¦   ¦   ¦--threads: bcftools_view_threads                                                                               
542  ¦   ¦   °--output_name: $(inputs.input1.basename.split(".filtered.snp.vcf")[0].concat(".concat.vcf"))                   
543  ¦   °--out:                                                                                                             
544  ¦       °--output                                                                                                       
545  ¦--bcftools_view_hard_filter                                                                                            
546  ¦   ¦--run:                                                                                                             
547  ¦   ¦   °--../wrappers/bcftools-view.cwl                                                                                
548  ¦   ¦--scatter:                                                                                                         
549  ¦   ¦   °--input                                                                                                        
550  ¦   ¦--in:                                                                                                              
551  ¦   ¦   ¦--input: bcftools_concat/output                                                                                
552  ¦   ¦   ¦--threads: bcftools_view_threads                                                                               
553  ¦   ¦   ¦--include: bcftools_view_include_hard_filters                                                                  
554  ¦   ¦   °--output_name: $(inputs.input.basename.split(".concat.vcf")[0].concat(".bcftools.hard.filtered.vcf"))          
555  ¦   °--out:                                                                                                             
556  ¦       °--output                                                                                                       
557  ¦--bcftools_norm_hard_filter                                                                                            
558  ¦   ¦--run:                                                                                                             
559  ¦   ¦   °--../wrappers/bcftools-norm.cwl                                                                                
560  ¦   ¦--scatter:                                                                                                         
561  ¦   ¦   °--input                                                                                                        
562  ¦   ¦--in:                                                                                                              
563  ¦   ¦   ¦--input: bcftools_view_hard_filter/output                                                                      
564  ¦   ¦   ¦--threads: bcftools_norm_threads                                                                               
565  ¦   ¦   ¦--reference: reference_genome                                                                                  
566  ¦   ¦   ¦--multiallelics: bcftoomls_norm_multiallelics                                                                  
567  ¦   ¦   °--output_type: $( "v" )                                                                                        
568  ¦   °--out:                                                                                                             
569  ¦       °--output                                                                                                       
570  ¦--table_annovar_hard_filtered                                                                                          
571  ¦   ¦--run:                                                                                                             
572  ¦   ¦   °--../wrappers/table-annovar.cwl                                                                                
573  ¦   ¦--scatter:                                                                                                         
574  ¦   ¦   °--query_file                                                                                                   
575  ¦   ¦--in:                                                                                                              
576  ¦   ¦   ¦--query_file: bcftools_norm_hard_filter/output                                                                 
577  ¦   ¦   ¦--database_location: table_annovar_database_location                                                           
578  ¦   ¦   ¦--build_over: table_annovar_build_over                                                                         
579  ¦   ¦   ¦--remove: table_annovar_remove                                                                                 
580  ¦   ¦   ¦--protocol: table_annovar_protocol                                                                             
581  ¦   ¦   ¦--operation: table_annovar_operation                                                                           
582  ¦   ¦   ¦--na_string: table_annovar_na_string                                                                           
583  ¦   ¦   ¦--vcfinput: table_annovar_vcfinput                                                                             
584  ¦   ¦   ¦--otherinfo: table_annovar_otherinfo                                                                           
585  ¦   ¦   °--convert_arg: table_annovar_convert_arg                                                                       
586  ¦   °--out:                                                                                                             
587  ¦       ¦--multianno_vcf                                                                                                
588  ¦       ¦--multianno_txt                                                                                                
589  ¦       °--avinput                                                                                                      
590  ¦--gatk_CNNScoreVariants                                                                                                
591  ¦   ¦--run:                                                                                                             
592  ¦   ¦   °--../wrappers/gatk-CNNScoreVariants.cwl                                                                        
593  ¦   ¦--scatter:                                                                                                         
594  ¦   ¦   ¦--variant                                                                                                      
595  ¦   ¦   °--aligned_reads                                                                                                
596  ¦   ¦--scatterMethod:                                                                                                   
597  ¦   ¦   °--dotproduct                                                                                                   
598  ¦   ¦--in:                                                                                                              
599  ¦   ¦   ¦--reference: reference_genome                                                                                  
600  ¦   ¦   ¦--variant: gatk_haplotypecaller_subworkflow/o_gatk_MergeVCFs                                                   
601  ¦   ¦   ¦--aligned_reads: samtools_index_2/alignments_with_index                                                        
602  ¦   ¦   ¦--OUTPUT: $(inputs.variant.basename.split(".vcf.gz")[0].concat(".cnn.vcf"))                                    
603  ¦   ¦   °--tensor_type: $( "read_tensor" )                                                                              
604  ¦   °--out:                                                                                                             
605  ¦       °--output                                                                                                       
606  ¦--gatk_FilterVariantTranches                                                                                           
607  ¦   ¦--run:                                                                                                             
608  ¦   ¦   °--../wrappers/gatk-FilterVariantTranches.cwl                                                                   
609  ¦   ¦--scatter:                                                                                                         
610  ¦   ¦   °--variant                                                                                                      
611  ¦   ¦--in:                                                                                                              
612  ¦   ¦   ¦--variant: gatk_CNNScoreVariants/output                                                                        
613  ¦   ¦   ¦--OUTPUT: $(inputs.variant.basename.split(".cnn.vcf")[0].concat(".cnn_filtered.vcf"))                          
614  ¦   ¦   ¦--resource_1: FilterVariantTranches_resource_1                                                                 
615  ¦   ¦   ¦--resource_2: FilterVariantTranches_resource_2                                                                 
616  ¦   ¦   °--resource_3: FilterVariantTranches_resource_3                                                                 
617  ¦   °--out:                                                                                                             
618  ¦       °--output                                                                                                       
619  ¦--bcftools_view_filter_cnn                                                                                             
620  ¦   ¦--run:                                                                                                             
621  ¦   ¦   °--../wrappers/bcftools-view.cwl                                                                                
622  ¦   ¦--scatter:                                                                                                         
623  ¦   ¦   °--input                                                                                                        
624  ¦   ¦--in:                                                                                                              
625  ¦   ¦   ¦--input: gatk_FilterVariantTranches/output                                                                     
626  ¦   ¦   ¦--threads: bcftools_view_threads                                                                               
627  ¦   ¦   ¦--include: bcftools_view_include_CNN_filters                                                                   
628  ¦   ¦   °--output_name: $(inputs.input.basename.split(".cnn_filtered.vcf")[0].concat(".bcftools_cnn_filtered.vcf"))     
629  ¦   °--out:                                                                                                             
630  ¦       °--output                                                                                                       
631  ¦--bcftools_norm_cnn                                                                                                    
632  ¦   ¦--run:                                                                                                             
633  ¦   ¦   °--../wrappers/bcftools-norm.cwl                                                                                
634  ¦   ¦--scatter:                                                                                                         
635  ¦   ¦   °--input                                                                                                        
636  ¦   ¦--in:                                                                                                              
637  ¦   ¦   ¦--input: bcftools_view_filter_cnn/output                                                                       
638  ¦   ¦   ¦--threads: bcftools_norm_threads                                                                               
639  ¦   ¦   ¦--reference: reference_genome                                                                                  
640  ¦   ¦   ¦--multiallelics: bcftoomls_norm_multiallelics                                                                  
641  ¦   ¦   °--output_type: $( "v" )                                                                                        
642  ¦   °--out:                                                                                                             
643  ¦       °--output                                                                                                       
644  °--table_annovar_cnn_filtered                                                                                           
645      ¦--run:                                                                                                             
646      ¦   °--../wrappers/table-annovar.cwl                                                                                
647      ¦--scatter:                                                                                                         
648      ¦   °--query_file                                                                                                   
649      ¦--in:                                                                                                              
650      ¦   ¦--query_file: bcftools_norm_cnn/output                                                                         
651      ¦   ¦--database_location: table_annovar_database_location                                                           
652      ¦   ¦--build_over: table_annovar_build_over                                                                         
653      ¦   ¦--remove: table_annovar_remove                                                                                 
654      ¦   ¦--protocol: table_annovar_protocol                                                                             
655      ¦   ¦--operation: table_annovar_operation                                                                           
656      ¦   ¦--na_string: table_annovar_na_string                                                                           
657      ¦   ¦--vcfinput: table_annovar_vcfinput                                                                             
658      ¦   ¦--otherinfo: table_annovar_otherinfo                                                                           
659      ¦   °--convert_arg: table_annovar_convert_arg                                                                       
660      °--out:                                                                                                             
661          ¦--multianno_vcf                                                                                                
662          ¦--multianno_txt                                                                                                
663          °--avinput                                                                                                      
```

</details>