# CWL-based ChIP-Seq workflow

## Description

A CWL-based pipeline for processing ChIP-Seq data (FASTQ format) and performing: 

- Peak calling
- Consensus peak count table generation
- Detection of super-enhancer regions
- Differential binding analysis

A pre-configured YAML template, based on validation analysis of publicly available HTS data, is available as example in the ``yaml_files`` folder. Moreover, tables of metadata, based on the same validation analysis, are available in the files ``EZH2_metadata_CLL.csv`` and ``H3K27me3_metadata_CLL.csv`` to serve as input examples for the design of comparisons during differential binding analysis. In addition, a list of ChIP-Seq blacklisted regions (human genome version 38; hg38) from the ENCODE project, which is can be used as input for the workflow, is provided in BED format (``hg38-blacklist.v2.bed``).

Briefly, the workflow performs the following steps:

1. Quality control of short reads (FastQC)
2. Trimming of the reads (e.g., removal of adapter and/or low quality sequences) (Trimmomatic)
3. Mapping to reference genome (HISAT2)
5. Convertion of mapped reads from SAM (Sequence Alignment Map) to BAM (Binary Alignment Map) format (samtools)
6. Sorting mapped reads based on chromosomal coordinates (samtools)
7. Adding information regarding paired end reads (e.g., CIGAR field information) (samtools)
8. Re-sorting based on chromosomal coordinates (samtools)
9. Removal of duplicate reads (samtools)
10. Index creation for coordinate-sorted BAM files to enable fast random access (samtools)
11. Production of quality metrics and files for the inspection of the mapped ChIP-Seq reads, taking into consideration the experimental design (deeptools2):
 - Read coverages for genomic regions of two or more BAM files are computed (multiBamSummary). The results are produced in compressed numpy array (NPZ) format and are used to calculate and visualize pairwise correlation values between the read coverages (plotCorrelation). 
 - Estimation of sequencing depth, through genomic position (base pair) sampling, and visualization is performed for multiple BAM files (plotCoverage).
 - Cumulative read coverages for each indexed BAM file are plotted by counting and sorting all reads overlapping a “window” of specified length (plotFingerprint).
 - Production of coverage track files (bigWig), with the coverage calculated as the number of reads per consecutive windows of predefined size (bamCoverage), and normalized through various available methods (e.g., Reads Per Kilobase per Million mapped reads; RPKM). The coverage track files are used to calculate scores per selected genomic regions (computeMatrix), typically genes, and a heatmap, based on the scores associated with these genomic regions, is produced (plotHeatmap).
12. Calling potential binding positions (peaks) to the genome (peak calling) (MACS2)
13. Generation of consensus peak count table for the application of custom analyses on MACS2 peak calling results (bedtools)
14. Detection of super-enhancer regions (Rank Ordering of Super-Enhancers; ROSE)
15. Differential binding analyses (DiffBind) for:
 - MACS2 peak calling results
 - ROSE-detected super-enhancer regions 

## CWL wrappers/software tools used in this pipeline

| Software | CWL.wrapper | CWL.type | Docker.image |
|---|---|---|---|
| - | get-raw-files.cwl | ExpressionTool | - |
| - | split-single-paired_v2.cwl | ExpressionTool | - |
| Trimmomatic | Trimmomatic_SE.cwl | CommandLineTool | staphb/trimmomatic:latest |
| Trimmomatic | Trimmomatic_PE.cwl | CommandLineTool | staphb/trimmomatic:latest |
| FastQC | fastqc.cwl | CommandLineTool | biowardrobe2/fastqc:v0.11.5 |
| cp | cp.cwl | CommandLineTool | ubuntu:latest |
| cp | rename.cwl | CommandLineTool | ubuntu:latest |
| - | check_trimming.cwl | ExpressionTool | - |
| HISAT2 | hisat2_v2.cwl | CommandLineTool | greatfireball/hisat2 |
| - | collect-hisat2-sam-files.cwl | ExpressionTool | - |
| SAMtools | samtools-view.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| SAMtools | samtools-sort.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| SAMtools | samtools-fixmate.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| SAMtools | samtools-markdup.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| SAMtools | samtools-index.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| Deeptools  | multiBamSummary.cwl | CommandLineTool | kerstenbreuer/deeptools:3.1.1 |
| Deeptools  | plotCorrelation.cwl | CommandLineTool | kerstenbreuer/deeptools:3.1.1 |
| Deeptools  | plotCoverage.cwl | CommandLineTool | kerstenbreuer/deeptools:3.1.1 |
| Deeptools  | plotFingerprint.cwl | CommandLineTool | kerstenbreuer/deeptools:3.1.1 |
| Deeptools  | bamCoverage.cwl | CommandLineTool | kerstenbreuer/deeptools:3.1.1 |
| Deeptools  | computeMatrix.cwl | CommandLineTool | kerstenbreuer/deeptools:3.1.1 |
| Deeptools  | plotHeatmap.cwl | CommandLineTool | kerstenbreuer/deeptools:3.1.1 |
| - | separate_control_treatment_files.cwl | ExpressionTool | - |
| MACS2 | macs2-callpeak.cwl | CommandLineTool | biowardrobe2/macs2:v2.1.1 |
| cat/awk/csvtk | table_total_peaks.cwl | CommandLineTool | quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0 |
| sort | sort.cwl | CommandLineTool | ubuntu:latest |
| bedtools | bedtools-merge.cwl | CommandLineTool | staphb/bedtools:2.30.0 |
| bedtools | bedtools-intersect.cwl | CommandLineTool | staphb/bedtools:2.30.0 |
| bedtools | bedtools-coverage.cwl | CommandLineTool | staphb/bedtools:2.30.0 |
| awk | awk.cwl | CommandLineTool | ubuntu:latest |
| paste | paste.cwl | CommandLineTool | ubuntu:latest |
| printf | printf.cwl | CommandLineTool | ubuntu:latest |
| ChIPQC | ChIPQC.cwl | CommandLineTool | biodataanalysisgroup/chipqc-diffbind-rscripts:v1.0 |
| DiffBind | DiffBind.cwl | CommandLineTool | biodataanalysisgroup/chipqc-diffbind-rscripts:v1.0 |
| bedtools | bedtools-intersect-narrowPeak.cwl | CommandLineTool | staphb/bedtools:2.30.0 |
| awk | bed_to_rose_gff.cwl | CommandLineTool | ubuntu:latest |
| ROSE | rose.cwl | CommandLineTool | biodataanalysisgroup/rose_main:v1.0 |
| awk | enhancer_bed_processing.cwl | CommandLineTool | ubuntu:latest |

## Workflow structure

A tree-based representation is available below for the inspection of all workflow steps. This tree-based hierarchical structure was produced using the [data.tree](https://cran.r-project.org/web/packages/data.tree/index.html) package.

```bash
1   CWL-based ChIP-Seq workflow                                                                                       
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
21   ¦       ¦--paired_files_fwd                                                                                      
22   ¦       ¦--paired_files_rev                                                                                      
23   ¦       ¦--fastqc_for_raw                                                                                        
24   ¦       ¦--fastqc_for_single                                                                                     
25   ¦       ¦--fastqc_for_paired                                                                                     
26   ¦       ¦--single_files_basenames                                                                                
27   ¦       ¦--paired_files_basenames                                                                                
28   ¦       ¦--single_files_sam                                                                                      
29   ¦       ¦--paired_files_sam                                                                                      
30   ¦       ¦--cp_command_raw                                                                                        
31   ¦       ¦--cp_command_single                                                                                     
32   ¦       ¦--cp_command_paired                                                                                     
33   ¦       ¦--trimmomatic_command_single                                                                            
34   ¦       °--trimmomatic_command_paired                                                                            
35   ¦--trimmomatic_single_end                                                                                        
36   ¦   ¦--run:                                                                                                      
37   ¦   ¦   °--../wrappers/Trimmomatic_SE.cwl                                                                        
38   ¦   ¦--scatter:                                                                                                  
39   ¦   ¦   °--input_fastq                                                                                           
40   ¦   ¦--in:                                                                                                       
41   ¦   ¦   ¦--input_fastq: split_single_paired/single_files                                                         
42   ¦   ¦   ¦--trimm_se_threads: trimmomatic_se_threads                                                              
43   ¦   ¦   ¦--command: split_single_paired/trimmomatic_command_single                                               
44   ¦   ¦   ¦--file_split: input_file_split                                                                          
45   ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                                    
46   ¦   ¦   ¦--illuminaClip: trimmomatic_se_illuminaClip                                                             
47   ¦   ¦   ¦--slidingWindow: trimmomatic_se_slidingWindow                                                           
48   ¦   ¦   ¦--leading: trimmomatic_se_leading                                                                       
49   ¦   ¦   ¦--trailing: trimmomatic_se_trailing                                                                     
50   ¦   ¦   °--minlen: trimmomatic_se_minlen                                                                         
51   ¦   °--out:                                                                                                      
52   ¦       ¦--stderr_log                                                                                            
53   ¦       °--outFastq                                                                                              
54   ¦--trimmomatic_paired_end                                                                                        
55   ¦   ¦--run:                                                                                                      
56   ¦   ¦   °--../wrappers/Trimmomatic_PE.cwl                                                                        
57   ¦   ¦--scatter:                                                                                                  
58   ¦   ¦   ¦--input_fastq_fwd                                                                                       
59   ¦   ¦   °--input_fastq_rev                                                                                       
60   ¦   ¦--scatterMethod:                                                                                            
61   ¦   ¦   °--dotproduct                                                                                            
62   ¦   ¦--in:                                                                                                       
63   ¦   ¦   ¦--input_fastq_fwd: split_single_paired/paired_files_fwd                                                 
64   ¦   ¦   ¦--input_fastq_rev: split_single_paired/paired_files_rev                                                 
65   ¦   ¦   ¦--trimm_pe_threads: trimmomatic_pe_threads                                                              
66   ¦   ¦   ¦--command: split_single_paired/trimmomatic_command_paired                                               
67   ¦   ¦   ¦--file_split: input_file_split                                                                          
68   ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                                    
69   ¦   ¦   ¦--file_split_rev: input_file_split_rev                                                                  
70   ¦   ¦   ¦--illuminaClip: trimmomatic_pe_illuminaClip                                                             
71   ¦   ¦   ¦--slidingWindow: trimmomatic_pe_slidingWindow                                                           
72   ¦   ¦   ¦--leading: trimmomatic_pe_leading                                                                       
73   ¦   ¦   ¦--trailing: trimmomatic_pe_trailing                                                                     
74   ¦   ¦   °--minlen: trimmomatic_pe_minlen                                                                         
75   ¦   °--out:                                                                                                      
76   ¦       ¦--stderr_log                                                                                            
77   ¦       ¦--outFastq_fwd_paired                                                                                   
78   ¦       ¦--outFastq_fwd_unpaired                                                                                 
79   ¦       ¦--outFastq_rev_paired                                                                                   
80   ¦       °--outFastq_rev_unpaired                                                                                 
81   ¦--fastqc_raw                                                                                                    
82   ¦   ¦--run:                                                                                                      
83   ¦   ¦   °--../wrappers/fastqc.cwl                                                                                
84   ¦   ¦--in:                                                                                                       
85   ¦   ¦   ¦--command: split_single_paired/fastqc_for_raw                                                           
86   ¦   ¦   °--input_files: get_raw_files/raw_files                                                                  
87   ¦   °--out:                                                                                                      
88   ¦       ¦--html_file                                                                                             
89   ¦       °--zipped_file                                                                                           
90   ¦--fastqc_single_trimmed                                                                                         
91   ¦   ¦--run:                                                                                                      
92   ¦   ¦   °--../wrappers/fastqc.cwl                                                                                
93   ¦   ¦--in:                                                                                                       
94   ¦   ¦   ¦--command: split_single_paired/fastqc_for_single                                                        
95   ¦   ¦   °--input_files: trimmomatic_single_end/outFastq                                                          
96   ¦   °--out:                                                                                                      
97   ¦       ¦--html_file                                                                                             
98   ¦       °--zipped_file                                                                                           
99   ¦--fastqc_paired_trimmed_fwd                                                                                     
100  ¦   ¦--run:                                                                                                      
101  ¦   ¦   °--../wrappers/fastqc.cwl                                                                                
102  ¦   ¦--in:                                                                                                       
103  ¦   ¦   ¦--command: split_single_paired/fastqc_for_paired                                                        
104  ¦   ¦   °--input_files: trimmomatic_paired_end/outFastq_fwd_paired                                               
105  ¦   °--out:                                                                                                      
106  ¦       ¦--html_file                                                                                             
107  ¦       °--zipped_file                                                                                           
108  ¦--fastqc_paired_trimmed_rev                                                                                     
109  ¦   ¦--run:                                                                                                      
110  ¦   ¦   °--../wrappers/fastqc.cwl                                                                                
111  ¦   ¦--in:                                                                                                       
112  ¦   ¦   ¦--command: split_single_paired/fastqc_for_paired                                                        
113  ¦   ¦   °--input_files: trimmomatic_paired_end/outFastq_rev_paired                                               
114  ¦   °--out:                                                                                                      
115  ¦       ¦--html_file                                                                                             
116  ¦       °--zipped_file                                                                                           
117  ¦--cp_fastqc_raw_zip                                                                                             
118  ¦   ¦--run:                                                                                                      
119  ¦   ¦   °--../wrappers/cp.cwl                                                                                    
120  ¦   ¦--in:                                                                                                       
121  ¦   ¦   ¦--command: split_single_paired/cp_command_raw                                                           
122  ¦   ¦   ¦--input_files: fastqc_raw/zipped_file                                                                   
123  ¦   ¦   °--outputdir: $( "fastqc_raw_zip" )                                                                      
124  ¦   °--out:                                                                                                      
125  ¦       °--output_dir                                                                                            
126  ¦--cp_fastqc_single_zip                                                                                          
127  ¦   ¦--run:                                                                                                      
128  ¦   ¦   °--../wrappers/cp.cwl                                                                                    
129  ¦   ¦--in:                                                                                                       
130  ¦   ¦   ¦--command: split_single_paired/cp_command_single                                                        
131  ¦   ¦   ¦--input_files: fastqc_single_trimmed/zipped_file                                                        
132  ¦   ¦   °--outputdir: $( "fastqc_single_trimmed_zip" )                                                           
133  ¦   °--out:                                                                                                      
134  ¦       °--output_dir                                                                                            
135  ¦--cp_fastqc_paired_zip                                                                                          
136  ¦   ¦--run:                                                                                                      
137  ¦   ¦   °--../wrappers/cp_paired.cwl                                                                             
138  ¦   ¦--in:                                                                                                       
139  ¦   ¦   ¦--command: split_single_paired/cp_command_paired                                                        
140  ¦   ¦   ¦--input_files_fwd: fastqc_paired_trimmed_fwd/zipped_file                                                
141  ¦   ¦   ¦--input_files_rev: fastqc_paired_trimmed_rev/zipped_file                                                
142  ¦   ¦   °--outputdir: $( "fastqc_paired_trimmed_zip" )                                                           
143  ¦   °--out:                                                                                                      
144  ¦       °--output_dir                                                                                            
145  ¦--rename_fastqc_raw_html                                                                                        
146  ¦   ¦--run:                                                                                                      
147  ¦   ¦   °--../wrappers/rename.cwl                                                                                
148  ¦   ¦--scatter:                                                                                                  
149  ¦   ¦   °--input_file                                                                                            
150  ¦   ¦--in:                                                                                                       
151  ¦   ¦   ¦--input_file: fastqc_raw/html_file                                                                      
152  ¦   ¦   °--run_type: $( "fastqc_raw_" )                                                                          
153  ¦   °--out:                                                                                                      
154  ¦       °--renamed_file                                                                                          
155  ¦--rename_fastqc_single_html                                                                                     
156  ¦   ¦--run:                                                                                                      
157  ¦   ¦   °--../wrappers/rename.cwl                                                                                
158  ¦   ¦--scatter:                                                                                                  
159  ¦   ¦   °--input_file                                                                                            
160  ¦   ¦--in:                                                                                                       
161  ¦   ¦   ¦--input_file: fastqc_single_trimmed/html_file                                                           
162  ¦   ¦   °--run_type: $( "fastqc_single_trimmed_" )                                                               
163  ¦   °--out:                                                                                                      
164  ¦       °--renamed_file                                                                                          
165  ¦--rename_fastqc_paired_html_fwd                                                                                 
166  ¦   ¦--run:                                                                                                      
167  ¦   ¦   °--../wrappers/rename.cwl                                                                                
168  ¦   ¦--scatter:                                                                                                  
169  ¦   ¦   °--input_file                                                                                            
170  ¦   ¦--in:                                                                                                       
171  ¦   ¦   ¦--input_file: fastqc_paired_trimmed_fwd/html_file                                                       
172  ¦   ¦   °--run_type: $( "fastqc_paired_trimmed_" )                                                               
173  ¦   °--out:                                                                                                      
174  ¦       °--renamed_file                                                                                          
175  ¦--rename_fastqc_paired_html_rev                                                                                 
176  ¦   ¦--run:                                                                                                      
177  ¦   ¦   °--../wrappers/rename.cwl                                                                                
178  ¦   ¦--scatter:                                                                                                  
179  ¦   ¦   °--input_file                                                                                            
180  ¦   ¦--in:                                                                                                       
181  ¦   ¦   ¦--input_file: fastqc_paired_trimmed_rev/html_file                                                       
182  ¦   ¦   °--run_type: $( "fastqc_paired_trimmed_" )                                                               
183  ¦   °--out:                                                                                                      
184  ¦       °--renamed_file                                                                                          
185  ¦--check_trimming                                                                                                
186  ¦   ¦--run:                                                                                                      
187  ¦   ¦   °--../wrappers/check_trimming.cwl                                                                        
188  ¦   ¦--in:                                                                                                       
189  ¦   ¦   ¦--trimming_check: input_trimming_check                                                                  
190  ¦   ¦   ¦--file_split: input_file_split                                                                          
191  ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                                    
192  ¦   ¦   ¦--file_split_rev: input_file_split_rev                                                                  
193  ¦   ¦   ¦--input_single: split_single_paired/single_files                                                        
194  ¦   ¦   ¦--input_paired_fwd: split_single_paired/paired_files_fwd                                                
195  ¦   ¦   ¦--input_paired_rev: split_single_paired/paired_files_rev                                                
196  ¦   ¦   ¦--trimming_single: trimmomatic_single_end/outFastq                                                      
197  ¦   ¦   ¦--trimming_paired_fwd: trimmomatic_paired_end/outFastq_fwd_paired                                       
198  ¦   ¦   °--trimming_paired_rev: trimmomatic_paired_end/outFastq_rev_paired                                       
199  ¦   °--out:                                                                                                      
200  ¦       ¦--single_files                                                                                          
201  ¦       ¦--paired_files_fwd                                                                                      
202  ¦       °--paired_files_rev                                                                                      
203  ¦--hisat2_for_single_reads                                                                                       
204  ¦   ¦--run:                                                                                                      
205  ¦   ¦   °--../wrappers/hisat2_v2.cwl                                                                             
206  ¦   ¦--scatter:                                                                                                  
207  ¦   ¦   ¦--files_with_unpaired_reads                                                                             
208  ¦   ¦   ¦--SAM_output                                                                                            
209  ¦   ¦   °--stderr_report                                                                                         
210  ¦   ¦--scatterMethod:                                                                                            
211  ¦   ¦   °--dotproduct                                                                                            
212  ¦   ¦--in:                                                                                                       
213  ¦   ¦   ¦--num_of_threads: hisat2_num_of_threads                                                                 
214  ¦   ¦   ¦--idx_directory: hisat2_idx_directory                                                                   
215  ¦   ¦   ¦--idx_basename: hisat2_idx_basename                                                                     
216  ¦   ¦   ¦--files_with_unpaired_reads: check_trimming/single_files                                                
217  ¦   ¦   ¦--SAM_output: split_single_paired/single_files_sam                                                      
218  ¦   ¦   °--stderr_report: split_single_paired/single_files_basenames                                             
219  ¦   °--out:                                                                                                      
220  ¦       ¦--output                                                                                                
221  ¦       °--output_stderr                                                                                         
222  ¦--hisat2_for_paired_reads                                                                                       
223  ¦   ¦--run:                                                                                                      
224  ¦   ¦   °--../wrappers/hisat2_v2.cwl                                                                             
225  ¦   ¦--scatter:                                                                                                  
226  ¦   ¦   ¦--files_with_first_mates                                                                                
227  ¦   ¦   ¦--files_with_second_mates                                                                               
228  ¦   ¦   ¦--SAM_output                                                                                            
229  ¦   ¦   °--stderr_report                                                                                         
230  ¦   ¦--scatterMethod:                                                                                            
231  ¦   ¦   °--dotproduct                                                                                            
232  ¦   ¦--in:                                                                                                       
233  ¦   ¦   ¦--num_of_threads: hisat2_num_of_threads                                                                 
234  ¦   ¦   ¦--idx_directory: hisat2_idx_directory                                                                   
235  ¦   ¦   ¦--idx_basename: hisat2_idx_basename                                                                     
236  ¦   ¦   ¦--files_with_first_mates: check_trimming/paired_files_fwd                                               
237  ¦   ¦   ¦--files_with_second_mates: check_trimming/paired_files_rev                                              
238  ¦   ¦   ¦--SAM_output: split_single_paired/paired_files_sam                                                      
239  ¦   ¦   °--stderr_report: split_single_paired/paired_files_basenames                                             
240  ¦   °--out:                                                                                                      
241  ¦       ¦--output                                                                                                
242  ¦       °--output_stderr                                                                                         
243  ¦--collect_hisat2_sam_files                                                                                      
244  ¦   ¦--run:                                                                                                      
245  ¦   ¦   °--../wrappers/collect-hisat2-sam-files.cwl                                                              
246  ¦   ¦--in:                                                                                                       
247  ¦   ¦   ¦--single_files: hisat2_for_single_reads/output                                                          
248  ¦   ¦   °--paired_files: hisat2_for_paired_reads/output                                                          
249  ¦   °--out:                                                                                                      
250  ¦       °--total_sam_files                                                                                       
251  ¦--samtools_view                                                                                                 
252  ¦   ¦--run:                                                                                                      
253  ¦   ¦   °--../wrappers/samtools-view.cwl                                                                         
254  ¦   ¦--scatter:                                                                                                  
255  ¦   ¦   °--input                                                                                                 
256  ¦   ¦--in:                                                                                                       
257  ¦   ¦   ¦--input: collect_hisat2_sam_files/total_sam_files                                                       
258  ¦   ¦   ¦--output_name: $( inputs.input.basename.split(".sam")[0].concat(".bam.tmp") )                           
259  ¦   ¦   ¦--readswithoutbits: samtools_readswithoutbits                                                           
260  ¦   ¦   ¦--samheader: $( true )                                                                                  
261  ¦   ¦   °--threads: samtools_view_threads                                                                        
262  ¦   °--out:                                                                                                      
263  ¦       °--output                                                                                                
264  ¦--samtools_sort_by_name                                                                                         
265  ¦   ¦--run:                                                                                                      
266  ¦   ¦   °--../wrappers/samtools-sort.cwl                                                                         
267  ¦   ¦--scatter:                                                                                                  
268  ¦   ¦   °--input                                                                                                 
269  ¦   ¦--in:                                                                                                       
270  ¦   ¦   ¦--compression_level: samtools_sort_compression_level                                                    
271  ¦   ¦   ¦--threads: samtools_sort_threads                                                                        
272  ¦   ¦   ¦--memory: samtools_sort_memory                                                                          
273  ¦   ¦   ¦--input: samtools_view/output                                                                           
274  ¦   ¦   ¦--output_name: $( inputs.input.basename.split(".bam.tmp")[0].concat(".name.sorted.bam") )               
275  ¦   ¦   °--sort_by_name: $( true )                                                                               
276  ¦   °--out:                                                                                                      
277  ¦       °--sorted                                                                                                
278  ¦--samtools_fixmate                                                                                              
279  ¦   ¦--run:                                                                                                      
280  ¦   ¦   °--../wrappers/samtools-fixmate.cwl                                                                      
281  ¦   ¦--scatter:                                                                                                  
282  ¦   ¦   °--input_file                                                                                            
283  ¦   ¦--in:                                                                                                       
284  ¦   ¦   ¦--threads: samtools_fixmate_threads                                                                     
285  ¦   ¦   ¦--output_format: samtools_fixmate_output_format                                                         
286  ¦   ¦   ¦--input_file: samtools_sort_by_name/sorted                                                              
287  ¦   ¦   °--output_file_name: $( inputs.input_file.basename.split(".name.sorted.bam")[0].concat("_unsorted.bam") )
288  ¦   °--out:                                                                                                      
289  ¦       °--output                                                                                                
290  ¦--samtools_sort                                                                                                 
291  ¦   ¦--run:                                                                                                      
292  ¦   ¦   °--../wrappers/samtools-sort.cwl                                                                         
293  ¦   ¦--scatter:                                                                                                  
294  ¦   ¦   °--input                                                                                                 
295  ¦   ¦--in:                                                                                                       
296  ¦   ¦   ¦--compression_level: samtools_sort_compression_level                                                    
297  ¦   ¦   ¦--threads: samtools_sort_threads                                                                        
298  ¦   ¦   ¦--memory: samtools_sort_memory                                                                          
299  ¦   ¦   ¦--input: samtools_fixmate/output                                                                        
300  ¦   ¦   °--output_name: $( inputs.input.basename.split("_unsorted.bam")[0].concat("_sorted.bam") )               
301  ¦   °--out:                                                                                                      
302  ¦       °--sorted                                                                                                
303  ¦--samtools_markdup                                                                                              
304  ¦   ¦--run:                                                                                                      
305  ¦   ¦   °--../wrappers/samtools-markdup.cwl                                                                      
306  ¦   ¦--scatter:                                                                                                  
307  ¦   ¦   °--input                                                                                                 
308  ¦   ¦--in:                                                                                                       
309  ¦   ¦   ¦--threads: samtools_markdup_threads                                                                     
310  ¦   ¦   ¦--remove_duplicates: $( true )                                                                          
311  ¦   ¦   ¦--input: samtools_sort/sorted                                                                           
312  ¦   ¦   °--output_name: $( inputs.input.basename.split("_sorted.bam")[0].concat("_markdup.bam") )                
313  ¦   °--out:                                                                                                      
314  ¦       °--output                                                                                                
315  ¦--samtools_index                                                                                                
316  ¦   ¦--run:                                                                                                      
317  ¦   ¦   °--../wrappers/samtools-index.cwl                                                                        
318  ¦   ¦--scatter:                                                                                                  
319  ¦   ¦   °--alignments                                                                                            
320  ¦   ¦--in:                                                                                                       
321  ¦   ¦   °--alignments: samtools_markdup/output                                                                   
322  ¦   °--out:                                                                                                      
323  ¦       °--alignments_with_index                                                                                 
324  ¦--multiBamSummary_file                                                                                          
325  ¦   ¦--run:                                                                                                      
326  ¦   ¦   °--../wrappers/multiBamSummary.cwl                                                                       
327  ¦   ¦--in:                                                                                                       
328  ¦   ¦   ¦--bam: samtools_index/alignments_with_index                                                             
329  ¦   ¦   ¦--threads: multiBamSummary_threads                                                                      
330  ¦   ¦   °--blackListFileName: blackListFile                                                                      
331  ¦   °--out:                                                                                                      
332  ¦       °--outNpz                                                                                                
333  ¦--plotCorrelation_file                                                                                          
334  ¦   ¦--run:                                                                                                      
335  ¦   ¦   °--../wrappers/plotCorrelation.cwl                                                                       
336  ¦   ¦--in:                                                                                                       
337  ¦   ¦   ¦--input: multiBamSummary_file/outNpz                                                                    
338  ¦   ¦   ¦--numbers: plotCorrelation_numbers                                                                      
339  ¦   ¦   ¦--method: plotCorrelation_method                                                                        
340  ¦   ¦   ¦--color: plotCorrelation_color                                                                          
341  ¦   ¦   ¦--title: plotCorrelation_title                                                                          
342  ¦   ¦   ¦--plotType: plotCorrelation_plotType                                                                    
343  ¦   ¦   °--outFileName: plotCorrelation_outFileName                                                              
344  ¦   °--out:                                                                                                      
345  ¦       °--outImage                                                                                              
346  ¦--plotCoverage_file                                                                                             
347  ¦   ¦--run:                                                                                                      
348  ¦   ¦   °--../wrappers/plotCoverage.cwl                                                                          
349  ¦   ¦--in:                                                                                                       
350  ¦   ¦   ¦--bam: samtools_index/alignments_with_index                                                             
351  ¦   ¦   ¦--threads: plotCoverage_threads                                                                         
352  ¦   ¦   ¦--skipZeros: $( true )                                                                                  
353  ¦   ¦   ¦--plotFileFormat: plotCoverage_plotFileFormat                                                           
354  ¦   ¦   ¦--outFileName: plotCoverage_outFileName                                                                 
355  ¦   ¦   °--blackListFileName: blackListFile                                                                      
356  ¦   °--out:                                                                                                      
357  ¦       °--outImage                                                                                              
358  ¦--plotFingerprint_file                                                                                          
359  ¦   ¦--run:                                                                                                      
360  ¦   ¦   °--../wrappers/plotFingerprint.cwl                                                                       
361  ¦   ¦--in:                                                                                                       
362  ¦   ¦   ¦--bam: samtools_index/alignments_with_index                                                             
363  ¦   ¦   ¦--threads: plotFingerprint_threads                                                                      
364  ¦   ¦   ¦--plotFileFormat: plotFingerprint_plotFileFormat                                                        
365  ¦   ¦   ¦--outFileName: plotFingerprint_outFileName                                                              
366  ¦   ¦   °--blackListFileName: blackListFile                                                                      
367  ¦   °--out:                                                                                                      
368  ¦       °--outImage                                                                                              
369  ¦--bamCoverage_norm                                                                                              
370  ¦   ¦--run:                                                                                                      
371  ¦   ¦   °--../wrappers/bamCoverage.cwl                                                                           
372  ¦   ¦--scatter:                                                                                                  
373  ¦   ¦   °--bam                                                                                                   
374  ¦   ¦--in:                                                                                                       
375  ¦   ¦   ¦--bam: samtools_index/alignments_with_index                                                             
376  ¦   ¦   ¦--effective_genome_size: bamCoverage_effective_genome_size                                              
377  ¦   ¦   ¦--normalizeUsing: bamCoverage_normalizeUsing                                                            
378  ¦   ¦   ¦--extendReads: bamCoverage_extendReads                                                                  
379  ¦   ¦   ¦--threads: bamCoverage_threads                                                                          
380  ¦   ¦   °--blackListFileName: blackListFile                                                                      
381  ¦   °--out:                                                                                                      
382  ¦       °--bigwig                                                                                                
383  ¦--computeMatrix                                                                                                 
384  ¦   ¦--run:                                                                                                      
385  ¦   ¦   °--../wrappers/computeMatrix.cwl                                                                         
386  ¦   ¦--in:                                                                                                       
387  ¦   ¦   ¦--bw: bamCoverage_norm/bigwig                                                                           
388  ¦   ¦   ¦--regions: computeMatrix_regions                                                                        
389  ¦   ¦   ¦--threads: computeMatrix_threads                                                                        
390  ¦   ¦   ¦--upstream: computeMatrix_upstream                                                                      
391  ¦   ¦   ¦--downstream: computeMatrix_downstream                                                                  
392  ¦   ¦   ¦--outputFile: computeMatrix_outputFile                                                                  
393  ¦   ¦   ¦--outFileSortedRegions: computeMatrix_outFileSortedRegions                                              
394  ¦   ¦   °--blackListFileName: blackListFile                                                                      
395  ¦   °--out:                                                                                                      
396  ¦       ¦--outMatrix                                                                                             
397  ¦       °--outRegions                                                                                            
398  ¦--plotHeatmap                                                                                                   
399  ¦   ¦--run:                                                                                                      
400  ¦   ¦   °--../wrappers/plotHeatmap.cwl                                                                           
401  ¦   ¦--in:                                                                                                       
402  ¦   ¦   ¦--matrix: computeMatrix/outMatrix                                                                       
403  ¦   ¦   ¦--plotFileFormat: plotHeatmap_plotFileFormat                                                            
404  ¦   ¦   °--outputFile: plotHeatmap_outputFile                                                                    
405  ¦   °--out:                                                                                                      
406  ¦       °--outHeatmap                                                                                            
407  ¦--separate_control_treatment_files                                                                              
408  ¦   ¦--run:                                                                                                      
409  ¦   ¦   °--../wrappers/separate_control_treatment_files.cwl                                                      
410  ¦   ¦--in:                                                                                                       
411  ¦   ¦   ¦--treatment_samples: input_treatment_samples                                                            
412  ¦   ¦   ¦--control_samples: input_control_samples                                                                
413  ¦   ¦   °--aligned_files: samtools_index/alignments_with_index                                                   
414  ¦   °--out:                                                                                                      
415  ¦       ¦--treatment_files                                                                                       
416  ¦       °--control_files                                                                                         
417  ¦--macs2_call_peaks                                                                                              
418  ¦   ¦--run:                                                                                                      
419  ¦   ¦   °--../wrappers/macs2-callpeak.cwl                                                                        
420  ¦   ¦--scatter:                                                                                                  
421  ¦   ¦   ¦--treatment                                                                                             
422  ¦   ¦   °--control                                                                                               
423  ¦   ¦--scatterMethod:                                                                                            
424  ¦   ¦   °--dotproduct                                                                                            
425  ¦   ¦--in:                                                                                                       
426  ¦   ¦   ¦--bdg: macs2_callpeak_bdg                                                                               
427  ¦   ¦   ¦--treatment: separate_control_treatment_files/treatment_files                                           
428  ¦   ¦   ¦--control: separate_control_treatment_files/control_files                                               
429  ¦   ¦   ¦--f: macs2_callpeak_format                                                                              
430  ¦   ¦   ¦--gsize: macs2_callpeak_gsize                                                                           
431  ¦   ¦   ¦--broad: macs2_callpeak_broad                                                                           
432  ¦   ¦   ¦--nomodel: macs2_callpeak_nomodel                                                                       
433  ¦   ¦   ¦--shift: macs2_shift                                                                                    
434  ¦   ¦   ¦--extsize: macs2_extsize                                                                                
435  ¦   ¦   ¦--p: macs2_pvalue                                                                                       
436  ¦   ¦   °--q: macs2_qvalue                                                                                       
437  ¦   °--out:                                                                                                      
438  ¦       ¦--narrowPeak                                                                                            
439  ¦       ¦--xls                                                                                                   
440  ¦       ¦--bed                                                                                                   
441  ¦       ¦--lambda                                                                                                
442  ¦       ¦--pileup                                                                                                
443  ¦       ¦--broadPeak                                                                                             
444  ¦       ¦--gappedPeak                                                                                            
445  ¦       ¦--model_r                                                                                               
446  ¦       °--cutoff                                                                                                
447  ¦--total_peaks_table                                                                                             
448  ¦   ¦--run:                                                                                                      
449  ¦   ¦   °--../wrappers/table_total_peaks.cwl                                                                     
450  ¦   ¦--in:                                                                                                       
451  ¦   ¦   ¦--input: macs2_call_peaks/narrowPeak                                                                    
452  ¦   ¦   °--consensus_bed: $( "narrowPeaks_consensus.bed" )                                                       
453  ¦   °--out:                                                                                                      
454  ¦       °--output                                                                                                
455  ¦--sort_peaks_table                                                                                              
456  ¦   ¦--run:                                                                                                      
457  ¦   ¦   °--../wrappers/sort.cwl                                                                                  
458  ¦   ¦--in:                                                                                                       
459  ¦   ¦   ¦--input: total_peaks_table/output                                                                       
460  ¦   ¦   ¦--key1: $("1,1")                                                                                        
461  ¦   ¦   ¦--key2: $("2,2n")                                                                                       
462  ¦   ¦   °--output_name: $( "narrowPeaks_consensus_sorted.bed" )                                                  
463  ¦   °--out:                                                                                                      
464  ¦       °--output                                                                                                
465  ¦--bedtools_merge                                                                                                
466  ¦   ¦--run:                                                                                                      
467  ¦   ¦   °--../wrappers/bedtools-merge.cwl                                                                        
468  ¦   ¦--in:                                                                                                       
469  ¦   ¦   ¦--input: sort_peaks_table/output                                                                        
470  ¦   ¦   °--output_name: $( "narrowPeaks_consensus_merged.bed" )                                                  
471  ¦   °--out:                                                                                                      
472  ¦       °--output                                                                                                
473  ¦--exclude_black_list_regions                                                                                    
474  ¦   ¦--run:                                                                                                      
475  ¦   ¦   °--../wrappers/bedtools-intersect.cwl                                                                    
476  ¦   ¦--in:                                                                                                       
477  ¦   ¦   ¦--feature_a: bedtools_merge/output                                                                      
478  ¦   ¦   ¦--feature_b: blackListFile                                                                              
479  ¦   ¦   ¦--no_overlap: $( true )                                                                                 
480  ¦   ¦   °--outputFile: $( "narrowPeaks_consensus_filtered.bed" )                                                 
481  ¦   °--out:                                                                                                      
482  ¦       °--output                                                                                                
483  ¦--bedtools_coverage                                                                                             
484  ¦   ¦--run:                                                                                                      
485  ¦   ¦   °--../wrappers/bedtools-coverage.cwl                                                                     
486  ¦   ¦--scatter:                                                                                                  
487  ¦   ¦   °--input_file                                                                                            
488  ¦   ¦--in:                                                                                                       
489  ¦   ¦   ¦--reference_file: exclude_black_list_regions/output                                                     
490  ¦   ¦   ¦--input_file: separate_control_treatment_files/treatment_files                                          
491  ¦   ¦   °--count_of_overlaps: $( true )                                                                          
492  ¦   °--out:                                                                                                      
493  ¦       °--output                                                                                                
494  ¦--extract_counts                                                                                                
495  ¦   ¦--run:                                                                                                      
496  ¦   ¦   °--../wrappers/awk.cwl                                                                                   
497  ¦   ¦--scatter:                                                                                                  
498  ¦   ¦   °--input_file                                                                                            
499  ¦   ¦--in:                                                                                                       
500  ¦   ¦   ¦--print: {print $4}                                                                                     
501  ¦   ¦   ¦--input_file: bedtools_coverage/output                                                                  
502  ¦   ¦   ¦--split_string: $("_markdup_counts.bed")                                                                
503  ¦   ¦   °--output_suffix: $(".count")                                                                            
504  ¦   °--out:                                                                                                      
505  ¦       °--output                                                                                                
506  ¦--extract_peaks                                                                                                 
507  ¦   ¦--run:                                                                                                      
508  ¦   ¦   °--../wrappers/awk.cwl                                                                                   
509  ¦   ¦--in:                                                                                                       
510  ¦   ¦   ¦--print: {print $1,$2,$3}                                                                               
511  ¦   ¦   ¦--input_file: exclude_black_list_regions/output                                                         
512  ¦   ¦   ¦--split_string: $(".bed")                                                                               
513  ¦   ¦   °--output_suffix: $(".coords")                                                                           
514  ¦   °--out:                                                                                                      
515  ¦       °--output                                                                                                
516  ¦--printf_header_samples                                                                                         
517  ¦   ¦--run:                                                                                                      
518  ¦   ¦   °--../wrappers/printf.cwl                                                                                
519  ¦   ¦--in:                                                                                                       
520  ¦   ¦   ¦--input_files: bedtools_coverage/output                                                                 
521  ¦   ¦   °--output_name: $( "final_consensus_count_samples.bed" )                                                 
522  ¦   °--out:                                                                                                      
523  ¦       °--output                                                                                                
524  ¦--paste_content_1                                                                                               
525  ¦   ¦--run:                                                                                                      
526  ¦   ¦   °--../wrappers/paste.cwl                                                                                 
527  ¦   ¦--in:                                                                                                       
528  ¦   ¦   ¦--input_files: extract_counts/output                                                                    
529  ¦   ¦   °--output_name: $( "consensus_count.txt" )                                                               
530  ¦   °--out:                                                                                                      
531  ¦       °--output                                                                                                
532  ¦--paste_content_2                                                                                               
533  ¦   ¦--run:                                                                                                      
534  ¦   ¦   °--../wrappers/paste.cwl                                                                                 
535  ¦   ¦--in:                                                                                                       
536  ¦   ¦   ¦--coordinates: extract_peaks/output                                                                     
537  ¦   ¦   ¦--counts: paste_content_1/output                                                                        
538  ¦   ¦   °--output_name: $( "consensus_count_with_coordinates.txt" )                                              
539  ¦   °--out:                                                                                                      
540  ¦       °--output                                                                                                
541  ¦--append_files                                                                                                  
542  ¦   ¦--run:                                                                                                      
543  ¦   ¦   °--../wrappers/awk.cwl                                                                                   
544  ¦   ¦--in:                                                                                                       
545  ¦   ¦   ¦--print: {print $0}                                                                                     
546  ¦   ¦   ¦--input_file: printf_header_samples/output                                                              
547  ¦   ¦   ¦--second_file: paste_content_2/output                                                                   
548  ¦   ¦   ¦--split_string: $("_samples.bed")                                                                       
549  ¦   ¦   °--output_suffix: $(".txt")                                                                              
550  ¦   °--out:                                                                                                      
551  ¦       °--output                                                                                                
552  ¦--ChIPQC_macs                                                                                                   
553  ¦   ¦--run:                                                                                                      
554  ¦   ¦   °--../wrappers/ChIPQC.cwl                                                                                
555  ¦   ¦--in:                                                                                                       
556  ¦   ¦   ¦--treatmentBAM: separate_control_treatment_files/treatment_files                                        
557  ¦   ¦   ¦--controlBAM: separate_control_treatment_files/control_files                                            
558  ¦   ¦   ¦--peaks: macs2_call_peaks/narrowPeak                                                                    
559  ¦   ¦   ¦--peakCaller: $("narrow")                                                                               
560  ¦   ¦   ¦--annotation: ChIPQC_annotation                                                                         
561  ¦   ¦   ¦--metadata: metadata_table                                                                              
562  ¦   ¦   ¦--consensus: ChIPQC_consensus                                                                           
563  ¦   ¦   ¦--bCount: ChIPQC_bCount                                                                                 
564  ¦   ¦   ¦--blacklist: ChIPQC_blacklist                                                                           
565  ¦   ¦   ¦--facetBy: ChIPQC_facetBy                                                                               
566  ¦   ¦   ¦--outputdir: $("ChIPQC_macs_HTML_report")                                                               
567  ¦   ¦   ¦--chipqc_experiment: $("ChIPQC_macs_experiment.rda")                                                    
568  ¦   ¦   °--chipqc_report: $("ChIPQC_macs_report.csv")                                                            
569  ¦   °--out:                                                                                                      
570  ¦       ¦--ChIPQCexperiment                                                                                      
571  ¦       ¦--outdir                                                                                                
572  ¦       °--ChIPQCreport                                                                                          
573  ¦--DiffBind_macs                                                                                                 
574  ¦   ¦--run:                                                                                                      
575  ¦   ¦   °--../wrappers/DiffBind.cwl                                                                              
576  ¦   ¦--in:                                                                                                       
577  ¦   ¦   ¦--treatmentBAM: separate_control_treatment_files/treatment_files                                        
578  ¦   ¦   ¦--controlBAM: separate_control_treatment_files/control_files                                            
579  ¦   ¦   ¦--peaks: macs2_call_peaks/narrowPeak                                                                    
580  ¦   ¦   ¦--peakCaller: $("narrow")                                                                               
581  ¦   ¦   ¦--metadata: metadata_table                                                                              
582  ¦   ¦   ¦--consensus: DiffBind_consensus                                                                         
583  ¦   ¦   ¦--minOverlap: DiffBind_minOverlap                                                                       
584  ¦   ¦   ¦--blacklist: DiffBind_blacklist                                                                         
585  ¦   ¦   ¦--greylist: DiffBind_greylist                                                                           
586  ¦   ¦   ¦--cores: DiffBind_cores                                                                                 
587  ¦   ¦   ¦--bParallel: DiffBind_bParallel                                                                         
588  ¦   ¦   ¦--normalization: DiffBind_normalization                                                                 
589  ¦   ¦   ¦--library: DiffBind_library                                                                             
590  ¦   ¦   ¦--background: DiffBind_background                                                                       
591  ¦   ¦   ¦--design: DiffBind_design                                                                               
592  ¦   ¦   ¦--reorderMeta_factor: DiffBind_reorderMeta_factor                                                       
593  ¦   ¦   ¦--reorderMeta_value: DiffBind_reorderMeta_value                                                         
594  ¦   ¦   ¦--retrieve_consensus: DiffBind_retrieve_consensus                                                       
595  ¦   ¦   ¦--low_read_count_filter: DiffBind_low_read_count_filter                                                 
596  ¦   ¦   ¦--diffbind_filterFun: DiffBind_filterFun                                                                
597  ¦   ¦   ¦--diffbind_results_filename: $("DiffBind_macs_diff_peaks.tsv")                                          
598  ¦   ¦   ¦--correlation_heatmap_filename: $("DiffBind_macs_diff_peaks_corr_heatmap.tif")                          
599  ¦   ¦   °--diffbind_consensus_filename: $("DiffBind_macs_diff_consensus_peaks.bed")                              
600  ¦   °--out:                                                                                                      
601  ¦       ¦--diffbind_results                                                                                      
602  ¦       ¦--correlation_heatmap                                                                                   
603  ¦       ¦--diffbind_consensus                                                                                    
604  ¦       ¦--diffbind_normalized_counts                                                                            
605  ¦       °--diffbind_dba_object                                                                                   
606  ¦--exclude_black_list_regions_narrowPeak                                                                         
607  ¦   ¦--run:                                                                                                      
608  ¦   ¦   °--../wrappers/bedtools-intersect-narrowPeak.cwl                                                         
609  ¦   ¦--scatter:                                                                                                  
610  ¦   ¦   °--feature_a                                                                                             
611  ¦   ¦--in:                                                                                                       
612  ¦   ¦   ¦--feature_a: macs2_call_peaks/narrowPeak                                                                
613  ¦   ¦   ¦--feature_b: blackListFile                                                                              
614  ¦   ¦   °--no_overlap: $( true )                                                                                 
615  ¦   °--out:                                                                                                      
616  ¦       °--output                                                                                                
617  ¦--bed_to_rose_gff_conversion                                                                                    
618  ¦   ¦--run:                                                                                                      
619  ¦   ¦   °--../wrappers/bed_to_rose_gff.cwl                                                                       
620  ¦   ¦--scatter:                                                                                                  
621  ¦   ¦   °--input                                                                                                 
622  ¦   ¦--in:                                                                                                       
623  ¦   ¦   °--input: exclude_black_list_regions_narrowPeak/output                                                   
624  ¦   °--out:                                                                                                      
625  ¦       °--output                                                                                                
626  ¦--rose_main                                                                                                     
627  ¦   ¦--run:                                                                                                      
628  ¦   ¦   °--../wrappers/rose.cwl                                                                                  
629  ¦   ¦--scatter:                                                                                                  
630  ¦   ¦   ¦--macs_peaks                                                                                            
631  ¦   ¦   ¦--ranking_bam                                                                                           
632  ¦   ¦   °--control_bam                                                                                           
633  ¦   ¦--scatterMethod:                                                                                            
634  ¦   ¦   °--dotproduct                                                                                            
635  ¦   ¦--in:                                                                                                       
636  ¦   ¦   ¦--macs_peaks: bed_to_rose_gff_conversion/output                                                         
637  ¦   ¦   ¦--ranking_bam: separate_control_treatment_files/treatment_files                                         
638  ¦   ¦   ¦--genome_build: rose_genome_build                                                                       
639  ¦   ¦   ¦--stitch_distance: rose_stitch_distance                                                                 
640  ¦   ¦   ¦--tss_distance: rose_tss_distance                                                                       
641  ¦   ¦   °--control_bam: separate_control_treatment_files/control_files                                           
642  ¦   °--out:                                                                                                      
643  ¦       ¦--gff_dir_outputs                                                                                       
644  ¦       ¦--mappedGFF_dir_outputs                                                                                 
645  ¦       ¦--STITCHED_ENHANCER_REGION_MAP                                                                          
646  ¦       ¦--AllEnhancers_table                                                                                    
647  ¦       ¦--SuperEnhancers_table                                                                                  
648  ¦       ¦--Enhancers_withSuper                                                                                   
649  ¦       °--Plot_points                                                                                           
650  ¦--enhancer_bed_processing                                                                                       
651  ¦   ¦--run:                                                                                                      
652  ¦   ¦   °--../wrappers/enhancer_bed_processing.cwl                                                               
653  ¦   ¦--scatter:                                                                                                  
654  ¦   ¦   °--input                                                                                                 
655  ¦   ¦--in:                                                                                                       
656  ¦   ¦   °--input: rose_main/Enhancers_withSuper                                                                  
657  ¦   °--out:                                                                                                      
658  ¦       °--output                                                                                                
659  ¦--ChIPQC_rose                                                                                                   
660  ¦   ¦--run:                                                                                                      
661  ¦   ¦   °--../wrappers/ChIPQC.cwl                                                                                
662  ¦   ¦--in:                                                                                                       
663  ¦   ¦   ¦--treatmentBAM: separate_control_treatment_files/treatment_files                                        
664  ¦   ¦   ¦--controlBAM: separate_control_treatment_files/control_files                                            
665  ¦   ¦   ¦--peaks: enhancer_bed_processing/output                                                                 
666  ¦   ¦   ¦--peakCaller: $("bed")                                                                                  
667  ¦   ¦   ¦--annotation: ChIPQC_annotation                                                                         
668  ¦   ¦   ¦--metadata: metadata_table                                                                              
669  ¦   ¦   ¦--consensus: ChIPQC_consensus                                                                           
670  ¦   ¦   ¦--bCount: ChIPQC_bCount                                                                                 
671  ¦   ¦   ¦--blacklist: ChIPQC_blacklist                                                                           
672  ¦   ¦   ¦--facetBy: ChIPQC_facetBy                                                                               
673  ¦   ¦   ¦--outputdir: $("ChIPQC_rose_HTML_report")                                                               
674  ¦   ¦   ¦--chipqc_experiment: $("ChIPQC_rose_experiment.rda")                                                    
675  ¦   ¦   °--chipqc_report: $("ChIPQC_rose_report.csv")                                                            
676  ¦   °--out:                                                                                                      
677  ¦       ¦--ChIPQCexperiment                                                                                      
678  ¦       ¦--outdir                                                                                                
679  ¦       °--ChIPQCreport                                                                                          
680  °--DiffBind_rose                                                                                                 
681      ¦--run:                                                                                                      
682      ¦   °--../wrappers/DiffBind.cwl                                                                              
683      ¦--in:                                                                                                       
684      ¦   ¦--treatmentBAM: separate_control_treatment_files/treatment_files                                        
685      ¦   ¦--controlBAM: separate_control_treatment_files/control_files                                            
686      ¦   ¦--peaks: enhancer_bed_processing/output                                                                 
687      ¦   ¦--peakCaller: $("bed")                                                                                  
688      ¦   ¦--metadata: metadata_table                                                                              
689      ¦   ¦--consensus: DiffBind_consensus                                                                         
690      ¦   ¦--minOverlap: DiffBind_minOverlap                                                                       
691      ¦   ¦--blacklist: DiffBind_blacklist                                                                         
692      ¦   ¦--greylist: DiffBind_greylist                                                                           
693      ¦   ¦--cores: DiffBind_cores                                                                                 
694      ¦   ¦--bParallel: DiffBind_bParallel                                                                         
695      ¦   ¦--normalization: DiffBind_normalization                                                                 
696      ¦   ¦--library: DiffBind_library                                                                             
697      ¦   ¦--background: DiffBind_background                                                                       
698      ¦   ¦--design: DiffBind_design                                                                               
699      ¦   ¦--reorderMeta_factor: DiffBind_reorderMeta_factor                                                       
700      ¦   ¦--reorderMeta_value: DiffBind_reorderMeta_value                                                         
701      ¦   ¦--retrieve_consensus: DiffBind_retrieve_consensus                                                       
702      ¦   ¦--low_read_count_filter: DiffBind_low_read_count_filter                                                 
703      ¦   ¦--diffbind_filterFun: DiffBind_filterFun                                                                
704      ¦   ¦--diffbind_results_filename: $("DiffBind_rose_diff_peaks.tsv")                                          
705      ¦   ¦--correlation_heatmap_filename: $("DiffBind_rose_diff_peaks_corr_heatmap.tif")                          
706      ¦   °--diffbind_consensus_filename: $("DiffBind_rose_diff_consensus_peaks.bed")                              
707      °--out:                                                                                                      
708          ¦--diffbind_results                                                                                      
709          ¦--correlation_heatmap                                                                                   
710          ¦--diffbind_consensus                                                                                    
711          ¦--diffbind_normalized_counts                                                                            
712          °--diffbind_dba_object                                                                                   
```
