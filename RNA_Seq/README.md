# CWL-based RNA-Seq workflow

## Description

A CWL-based pipeline for processing RNA-Seq data (FASTQ format) and performing differential gene/transcript expression analysis. 

Briefly, the workflow performs the following steps:

1. Quality control of Illumina reads (FastQC)
2. trimming of the reads (e.g., removal of adapter and/or low quality sequences) (Trim galore)
3. (Optional)  custom processing of the reads using FASTA/Q Trimmer (part of the FASTX-toolkit) 
4. Mapping to reference genome (HISAT2)
5. Convertion of mapped reads from SAM (Sequence Alignment Map) to BAM (Binary Alignment Map) format
6. Sorting mapped reads based on chromosomal coordinates

Subsequently, two independent workflows are implemented for differential expression analysis at the transcript and gene level. 

**First**, following the [reference protocol](https://doi.org/10.1038/nprot.2016.095) for HISAT, StringTie and Ballgown transcript expression analysis, StringTie along with a reference transcript annotation GTF (Gene Transfer Format) file (if one is available) is used to:

- Assemble transcripts for each RNA-Seq sample using the previous read alignments (BAM files)
- Generate a global, non-redundant set of transcripts observed in any of the RNA-Seq samples
- Estimate transcript abundances and generate read coverage tables for each RNA-Seq sample, based on the global, merged set of transcripts (rather than the reference) which is observed across all samples

Ballgown program is then used to load the coverage tables generated in the previous step and perform statistical analyses for differential expression at the transcript level. Notably, the StringTie - Ballgown protocol applied here was selected to include potentially novel transcripts in the analysis. 

**Second**, featureCounts is used to count reads that are mapped to selected genomic features, in this case genes by default, and generate a table of read counts per gene and sample. This table is passed as input to DESeq2 to perform differential expression analysis at the gene level. Both Ballgown and DESeq2 R scripts, along with their respective CWL wrappers, were designed to receive as input various parameters, such as experimental design, contrasts of interest, numeric thresholds, and hidden batch effects.

## CWL wrappers/software tools used in this pipeline

| Software | CWL.wrapper | CWL.type | Docker.image |
|---|---|---|---|
| - | get-raw-files.cwl | ExpressionTool | - |
| - | split-single-paired.cwl | ExpressionTool | - |
| Trim galore | trim-galore.cwl | CommandLineTool | kerstenbreuer/trim_galore:0.4.4_1.14_0.11.7 |
| FastQC | fastqc.cwl | CommandLineTool | biowardrobe2/fastqc:v0.11.5 |
| cp | cp.cwl | CommandLineTool | ubuntu:latest |
| cp | rename.cwl | CommandLineTool | ubuntu:latest |
| fastx_trimmer | fastx-trimmer.cwl | CommandLineTool | biowardrobe2/fastx_toolkit:v0.0.14 |
| - | check-selected-steps.cwl | ExpressionTool | - |
| HISAT2 | hisat2.cwl | CommandLineTool | greatfireball/hisat2 |
| - | collect-hisat2-sam-files.cwl | ExpressionTool | - |
| SAMtools | samtools-view.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| SAMtools | samtools-sort.cwl | CommandLineTool | quay.io/biocontainers/samtools:1.14--hb421002_0 |
| Stringtie | stringtie.cwl | CommandLineTool | gawbul/docker-stringtie |
| Stringtie | stringtie-for-ballgown.cwl | CommandLineTool | gawbul/docker-stringtie |
| Ballgown | ballgown.cwl | CommandLineTool | biodataanalysisgroup/ballgown-rscript:v1.0 |
| featureCounts | featureCounts.cwl | CommandLineTool | genomicpariscentre/featurecounts |
| DESeq2 | DESeq2.cwl | CommandLineTool | biodataanalysisgroup/deseq2-rscript:v1.0 |

## Workflow structure

A tree-based representation is available below for the inspection of all workflow steps. This tree-based hierarchical structure was produced using the [data.tree](https://cran.r-project.org/web/packages/data.tree/index.html) package.

```bash
1   CWL-based RNA-Seq workflow                                                                         
2    ¦--get_raw_files                                                                                  
3    ¦   ¦--run:                                                                                       
4    ¦   ¦   °--../wrappers/get-raw-files.cwl                                                          
5    ¦   ¦--in:                                                                                        
6    ¦   ¦   °--DIRECTORY: raw_files_directory                                                         
7    ¦   °--out:                                                                                       
8    ¦       °--raw_files                                                                              
9    ¦--split_single_paired                                                                            
10   ¦   ¦--run:                                                                                       
11   ¦   ¦   °--../wrappers/split-single-paired.cwl                                                    
12   ¦   ¦--in:                                                                                        
13   ¦   ¦   ¦--input_files: get_raw_files/raw_files                                                   
14   ¦   ¦   ¦--file_split: input_file_split                                                           
15   ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                     
16   ¦   ¦   ¦--file_split_rev: input_file_split_rev                                                   
17   ¦   ¦   ¦--qc_check: input_qc_check                                                               
18   ¦   ¦   °--trimming_check: input_trimming_check                                                   
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
29   ¦       ¦--cp_command_paired                                                                      
30   ¦       ¦--fastx_command_single                                                                   
31   ¦       °--fastx_command_paired                                                                   
32   ¦--trim_galore_single                                                                             
33   ¦   ¦--run:                                                                                       
34   ¦   ¦   °--../wrappers/trim-galore.cwl                                                            
35   ¦   ¦--in:                                                                                        
36   ¦   ¦   ¦--command: split_single_paired/trim_galore_for_single                                    
37   ¦   ¦   ¦--fq_files: split_single_paired/single_files                                             
38   ¦   ¦   ¦--length: tg_length                                                                      
39   ¦   ¦   ¦--quality: tg_quality                                                                    
40   ¦   ¦   ¦--compression: tg_compression                                                            
41   ¦   ¦   ¦--do_not_compress: tg_do_not_compress                                                    
42   ¦   ¦   ¦--trim_suffix: tg_trim_suffix                                                            
43   ¦   ¦   ¦--strigency: tg_strigency                                                                
44   ¦   ¦   °--paired: $( false )                                                                     
45   ¦   °--out:                                                                                       
46   ¦       ¦--trim_galore                                                                            
47   ¦       °--trim_galore_report                                                                     
48   ¦--trim_galore_paired                                                                             
49   ¦   ¦--run:                                                                                       
50   ¦   ¦   °--../wrappers/trim-galore.cwl                                                            
51   ¦   ¦--in:                                                                                        
52   ¦   ¦   ¦--command: split_single_paired/trim_galore_for_paired                                    
53   ¦   ¦   ¦--fq_files: split_single_paired/paired_files                                             
54   ¦   ¦   ¦--length: tg_length                                                                      
55   ¦   ¦   ¦--quality: tg_quality                                                                    
56   ¦   ¦   ¦--compression: tg_compression                                                            
57   ¦   ¦   ¦--do_not_compress: tg_do_not_compress                                                    
58   ¦   ¦   ¦--trim_suffix: tg_trim_suffix                                                            
59   ¦   ¦   ¦--strigency: tg_strigency                                                                
60   ¦   ¦   °--paired: $( true )                                                                      
61   ¦   °--out:                                                                                       
62   ¦       ¦--trim_galore                                                                            
63   ¦       °--trim_galore_report                                                                     
64   ¦--fastqc_raw                                                                                     
65   ¦   ¦--run:                                                                                       
66   ¦   ¦   °--../wrappers/fastqc.cwl                                                                 
67   ¦   ¦--in:                                                                                        
68   ¦   ¦   ¦--command: split_single_paired/fastqc_for_raw                                            
69   ¦   ¦   °--input_files: get_raw_files/raw_files                                                   
70   ¦   °--out:                                                                                       
71   ¦       ¦--html_file                                                                              
72   ¦       °--zipped_file                                                                            
73   ¦--fastqc_single_trimmed                                                                          
74   ¦   ¦--run:                                                                                       
75   ¦   ¦   °--../wrappers/fastqc.cwl                                                                 
76   ¦   ¦--in:                                                                                        
77   ¦   ¦   ¦--command: split_single_paired/fastqc_for_single                                         
78   ¦   ¦   °--input_files: trim_galore_single/trim_galore                                            
79   ¦   °--out:                                                                                       
80   ¦       ¦--html_file                                                                              
81   ¦       °--zipped_file                                                                            
82   ¦--fastqc_paired_trimmed                                                                          
83   ¦   ¦--run:                                                                                       
84   ¦   ¦   °--../wrappers/fastqc.cwl                                                                 
85   ¦   ¦--in:                                                                                        
86   ¦   ¦   ¦--command: split_single_paired/fastqc_for_paired                                         
87   ¦   ¦   °--input_files: trim_galore_paired/trim_galore                                            
88   ¦   °--out:                                                                                       
89   ¦       ¦--html_file                                                                              
90   ¦       °--zipped_file                                                                            
91   ¦--cp_fastqc_raw_zip                                                                              
92   ¦   ¦--run:                                                                                       
93   ¦   ¦   °--../wrappers/cp.cwl                                                                     
94   ¦   ¦--in:                                                                                        
95   ¦   ¦   ¦--command: split_single_paired/cp_command_raw                                            
96   ¦   ¦   ¦--input_files: fastqc_raw/zipped_file                                                    
97   ¦   ¦   °--outputdir: $( "fastqc_raw_zip" )                                                       
98   ¦   °--out:                                                                                       
99   ¦       °--output_dir                                                                             
100  ¦--cp_fastqc_single_zip                                                                           
101  ¦   ¦--run:                                                                                       
102  ¦   ¦   °--../wrappers/cp.cwl                                                                     
103  ¦   ¦--in:                                                                                        
104  ¦   ¦   ¦--command: split_single_paired/cp_command_single                                         
105  ¦   ¦   ¦--input_files: fastqc_single_trimmed/zipped_file                                         
106  ¦   ¦   °--outputdir: $( "fastqc_single_trimmed_zip" )                                            
107  ¦   °--out:                                                                                       
108  ¦       °--output_dir                                                                             
109  ¦--cp_fastqc_paired_zip                                                                           
110  ¦   ¦--run:                                                                                       
111  ¦   ¦   °--../wrappers/cp.cwl                                                                     
112  ¦   ¦--in:                                                                                        
113  ¦   ¦   ¦--command: split_single_paired/cp_command_paired                                         
114  ¦   ¦   ¦--input_files: fastqc_paired_trimmed/zipped_file                                         
115  ¦   ¦   °--outputdir: $( "fastqc_paired_trimmed_zip" )                                            
116  ¦   °--out:                                                                                       
117  ¦       °--output_dir                                                                             
118  ¦--rename_fastqc_raw_html                                                                         
119  ¦   ¦--run:                                                                                       
120  ¦   ¦   °--../wrappers/rename.cwl                                                                 
121  ¦   ¦--scatter:                                                                                   
122  ¦   ¦   °--input_file                                                                             
123  ¦   ¦--in:                                                                                        
124  ¦   ¦   ¦--input_file: fastqc_raw/html_file                                                       
125  ¦   ¦   °--run_type: $( "fastqc_raw_" )                                                           
126  ¦   °--out:                                                                                       
127  ¦       °--renamed_file                                                                           
128  ¦--rename_fastqc_single_html                                                                      
129  ¦   ¦--run:                                                                                       
130  ¦   ¦   °--../wrappers/rename.cwl                                                                 
131  ¦   ¦--scatter:                                                                                   
132  ¦   ¦   °--input_file                                                                             
133  ¦   ¦--in:                                                                                        
134  ¦   ¦   ¦--input_file: fastqc_single_trimmed/html_file                                            
135  ¦   ¦   °--run_type: $( "fastqc_single_trimmed_" )                                                
136  ¦   °--out:                                                                                       
137  ¦       °--renamed_file                                                                           
138  ¦--rename_fastqc_paired_html                                                                      
139  ¦   ¦--run:                                                                                       
140  ¦   ¦   °--../wrappers/rename.cwl                                                                 
141  ¦   ¦--scatter:                                                                                   
142  ¦   ¦   °--input_file                                                                             
143  ¦   ¦--in:                                                                                        
144  ¦   ¦   ¦--input_file: fastqc_paired_trimmed/html_file                                            
145  ¦   ¦   °--run_type: $( "fastqc_paired_trimmed_" )                                                
146  ¦   °--out:                                                                                       
147  ¦       °--renamed_file                                                                           
148  ¦--fastx_trimmer_single                                                                           
149  ¦   ¦--run:                                                                                       
150  ¦   ¦   °--../wrappers/fastx-trimmer.cwl                                                          
151  ¦   ¦--scatter:                                                                                   
152  ¦   ¦   °--input_file                                                                             
153  ¦   ¦--in:                                                                                        
154  ¦   ¦   ¦--command: split_single_paired/fastx_command_single                                      
155  ¦   ¦   ¦--first_base_to_keep: fastx_first_base_to_keep                                           
156  ¦   ¦   ¦--last_base_to_keep: fastx_last_base_to_keep                                             
157  ¦   ¦   °--input_file: trim_galore_single/trim_galore                                             
158  ¦   °--out:                                                                                       
159  ¦       °--output                                                                                 
160  ¦--fastx_trimmer_paired                                                                           
161  ¦   ¦--run:                                                                                       
162  ¦   ¦   °--../wrappers/fastx-trimmer.cwl                                                          
163  ¦   ¦--scatter:                                                                                   
164  ¦   ¦   °--input_file                                                                             
165  ¦   ¦--in:                                                                                        
166  ¦   ¦   ¦--command: split_single_paired/fastx_command_paired                                      
167  ¦   ¦   ¦--first_base_to_keep: fastx_first_base_to_keep                                           
168  ¦   ¦   ¦--last_base_to_keep: fastx_last_base_to_keep                                             
169  ¦   ¦   °--input_file: trim_galore_paired/trim_galore                                             
170  ¦   °--out:                                                                                       
171  ¦       °--output                                                                                 
172  ¦--check_for_fastx_and_produce_names                                                              
173  ¦   ¦--run:                                                                                       
174  ¦   ¦   °--../wrappers/check-selected-steps.cwl                                                   
175  ¦   ¦--in:                                                                                        
176  ¦   ¦   ¦--input_check: premapping_input_check                                                    
177  ¦   ¦   ¦--single_files: split_single_paired/single_files                                         
178  ¦   ¦   ¦--paired_files: split_single_paired/paired_files                                         
179  ¦   ¦   ¦--trim_galore_single: trim_galore_single/trim_galore                                     
180  ¦   ¦   ¦--trim_galore_paired: trim_galore_paired/trim_galore                                     
181  ¦   ¦   ¦--fastx_trimmer_single: fastx_trimmer_single/output                                      
182  ¦   ¦   ¦--fastx_trimmer_paired: fastx_trimmer_paired/output                                      
183  ¦   ¦   ¦--file_split: input_file_split                                                           
184  ¦   ¦   ¦--file_split_fwd_single: input_file_split_fwd_single                                     
185  ¦   ¦   °--file_split_rev: input_file_split_rev                                                   
186  ¦   °--out:                                                                                       
187  ¦       ¦--single_trim                                                                            
188  ¦       ¦--single_hisat2_sam                                                                      
189  ¦       ¦--paired_trim_1                                                                          
190  ¦       ¦--paired_trim_2                                                                          
191  ¦       °--paired_hisat2_sam                                                                      
192  ¦--hisat2_for_single_reads                                                                        
193  ¦   ¦--run:                                                                                       
194  ¦   ¦   °--../wrappers/hisat2.cwl                                                                 
195  ¦   ¦--scatter:                                                                                   
196  ¦   ¦   ¦--files_with_unpaired_reads                                                              
197  ¦   ¦   °--SAM_output                                                                             
198  ¦   ¦--scatterMethod:                                                                             
199  ¦   ¦   °--dotproduct                                                                             
200  ¦   ¦--in:                                                                                        
201  ¦   ¦   ¦--num_of_threads: hisat2_num_of_threads                                                  
202  ¦   ¦   ¦--alignments_tailored_trans_assemb: hisat2_alignments_tailored_trans_assemb              
203  ¦   ¦   ¦--idx_directory: hisat2_idx_directory                                                    
204  ¦   ¦   ¦--idx_basename: hisat2_idx_basename                                                      
205  ¦   ¦   ¦--files_with_unpaired_reads: check_for_fastx_and_produce_names/single_trim               
206  ¦   ¦   °--SAM_output: check_for_fastx_and_produce_names/single_hisat2_sam                        
207  ¦   °--out:                                                                                       
208  ¦       ¦--output                                                                                 
209  ¦       °--output_stderr                                                                          
210  ¦--hisat2_for_paired_reads                                                                        
211  ¦   ¦--run:                                                                                       
212  ¦   ¦   °--../wrappers/hisat2.cwl                                                                 
213  ¦   ¦--scatter:                                                                                   
214  ¦   ¦   ¦--files_with_first_mates                                                                 
215  ¦   ¦   ¦--files_with_second_mates                                                                
216  ¦   ¦   °--SAM_output                                                                             
217  ¦   ¦--scatterMethod:                                                                             
218  ¦   ¦   °--dotproduct                                                                             
219  ¦   ¦--in:                                                                                        
220  ¦   ¦   ¦--num_of_threads: hisat2_num_of_threads                                                  
221  ¦   ¦   ¦--known_splicesite_infile: hisat2_known_splicesite_infile                                
222  ¦   ¦   ¦--alignments_tailored_trans_assemb: hisat2_alignments_tailored_trans_assemb              
223  ¦   ¦   ¦--idx_directory: hisat2_idx_directory                                                    
224  ¦   ¦   ¦--idx_basename: hisat2_idx_basename                                                      
225  ¦   ¦   ¦--files_with_first_mates: check_for_fastx_and_produce_names/paired_trim_1                
226  ¦   ¦   ¦--files_with_second_mates: check_for_fastx_and_produce_names/paired_trim_2               
227  ¦   ¦   °--SAM_output: check_for_fastx_and_produce_names/paired_hisat2_sam                        
228  ¦   °--out:                                                                                       
229  ¦       ¦--output                                                                                 
230  ¦       °--output_stderr                                                                          
231  ¦--collect_hisat2_sam_files                                                                       
232  ¦   ¦--run:                                                                                       
233  ¦   ¦   °--../wrappers/collect-hisat2-sam-files.cwl                                               
234  ¦   ¦--in:                                                                                        
235  ¦   ¦   ¦--single_files: hisat2_for_single_reads/output                                           
236  ¦   ¦   °--paired_files: hisat2_for_paired_reads/output                                           
237  ¦   °--out:                                                                                       
238  ¦       °--total_sam_files                                                                        
239  ¦--samtools_view                                                                                  
240  ¦   ¦--run:                                                                                       
241  ¦   ¦   °--../wrappers/samtools-view.cwl                                                          
242  ¦   ¦--scatter:                                                                                   
243  ¦   ¦   °--input                                                                                  
244  ¦   ¦--in:                                                                                        
245  ¦   ¦   ¦--input: collect_hisat2_sam_files/total_sam_files                                        
246  ¦   ¦   ¦--output_name: $( inputs.input.basename.split(".sam")[0].concat("_unsorted.bam") )       
247  ¦   ¦   ¦--threads: samtools_view_threads                                                         
248  ¦   ¦   ¦--isbam: samtools_view_isbam                                                             
249  ¦   ¦   ¦--collapsecigar: samtools_view_collapsecigar                                             
250  ¦   ¦   ¦--uncompressed: samtools_view_uncompressed                                               
251  ¦   ¦   ¦--fastcompression: samtools_view_fastcompression                                         
252  ¦   ¦   ¦--samheader: samtools_view_samheader                                                     
253  ¦   ¦   ¦--count: samtools_view_count                                                             
254  ¦   ¦   ¦--readswithoutbits: samtools_view_readswithoutbits                                       
255  ¦   ¦   ¦--readsingroup: samtools_view_readsingroup                                               
256  ¦   ¦   ¦--readtagtostrip: samtools_view_readtagtostrip                                           
257  ¦   ¦   ¦--readsquality: samtools_view_readsquality                                               
258  ¦   ¦   ¦--readswithbits: samtools_view_readswithbits                                             
259  ¦   ¦   ¦--cigar: samtools_view_cigar                                                             
260  ¦   ¦   ¦--iscram: samtools_view_iscram                                                           
261  ¦   ¦   ¦--randomseed: samtools_view_randomseed                                                   
262  ¦   ¦   ¦--region: samtools_view_region                                                           
263  ¦   ¦   °--readsinlibrary: samtools_view_readsinlibrary                                           
264  ¦   °--out:                                                                                       
265  ¦       °--output                                                                                 
266  ¦--samtools_sort                                                                                  
267  ¦   ¦--run:                                                                                       
268  ¦   ¦   °--../wrappers/samtools-sort.cwl                                                          
269  ¦   ¦--scatter:                                                                                   
270  ¦   ¦   °--input                                                                                  
271  ¦   ¦--in:                                                                                        
272  ¦   ¦   ¦--input: samtools_view/output                                                            
273  ¦   ¦   ¦--compression_level: samtools_sort_compression_level                                     
274  ¦   ¦   ¦--threads: samtools_sort_threads                                                         
275  ¦   ¦   ¦--memory: samtools_sort_memory                                                           
276  ¦   ¦   ¦--output_name: $( inputs.input.basename.split("_unsorted.bam")[0].concat("_sorted.bam") )
277  ¦   ¦   °--sort_by_name: samtools_sort_sort_by_name                                               
278  ¦   °--out:                                                                                       
279  ¦       °--sorted                                                                                 
280  ¦--stringtie_transcript_assembly                                                                  
281  ¦   ¦--run:                                                                                       
282  ¦   ¦   °--../wrappers/stringtie.cwl                                                              
283  ¦   ¦--scatter:                                                                                   
284  ¦   ¦   °--input_bam                                                                              
285  ¦   ¦--in:                                                                                        
286  ¦   ¦   ¦--input_bam: samtools_sort/sorted                                                        
287  ¦   ¦   ¦--guide_gff: stringtie_guide_gff                                                         
288  ¦   ¦   ¦--sample_label: $( inputs.input_bam.basename.split("_sorted.bam")[0] )                   
289  ¦   ¦   ¦--cpus: stringtie_cpus                                                                   
290  ¦   ¦   ¦--out_gtf: $( inputs.input_bam.basename.split("_sorted.bam")[0].concat(".gtf") )         
291  ¦   ¦   ¦--verbose: stringtie_verbose                                                             
292  ¦   ¦   ¦--min_isoform_abundance: stringtie_min_isoform_abundance                                 
293  ¦   ¦   ¦--junction_coverage: stringtie_junction_coverage                                         
294  ¦   ¦   ¦--min_read_coverage: stringtie_min_read_coverage                                         
295  ¦   ¦   °--conservative_mode: stringtie_conservative_mode                                         
296  ¦   °--out:                                                                                       
297  ¦       °--output_gtf                                                                             
298  ¦--stringtie_merge                                                                                
299  ¦   ¦--run:                                                                                       
300  ¦   ¦   °--../wrappers/stringtie.cwl                                                              
301  ¦   ¦--scatterMethod:                                                                             
302  ¦   ¦   °--dotproduct                                                                             
303  ¦   ¦--in:                                                                                        
304  ¦   ¦   ¦--input_gtfs: stringtie_transcript_assembly/output_gtf                                   
305  ¦   ¦   ¦--transcript_merge_mode: stringtie_transcript_merge_mode                                 
306  ¦   ¦   ¦--guide_gff: stringtie_guide_gff                                                         
307  ¦   ¦   ¦--cpus: stringtie_cpus                                                                   
308  ¦   ¦   ¦--verbose: stringtie_verbose                                                             
309  ¦   ¦   °--out_gtf: stringtie_out_gtf                                                             
310  ¦   °--out:                                                                                       
311  ¦       °--output_gtf                                                                             
312  ¦--stringtie_expression                                                                           
313  ¦   ¦--run:                                                                                       
314  ¦   ¦   °--../wrappers/stringtie-for-ballgown.cwl                                                 
315  ¦   ¦--scatter:                                                                                   
316  ¦   ¦   °--input_bam                                                                              
317  ¦   ¦--in:                                                                                        
318  ¦   ¦   ¦--input_bam: samtools_sort/sorted                                                        
319  ¦   ¦   ¦--guide_gff: stringtie_merge/output_gtf                                                  
320  ¦   ¦   ¦--sample_label: $( inputs.input_bam.basename.split("_sorted.bam")[0] )                   
321  ¦   ¦   ¦--cpus: stringtie_cpus                                                                   
322  ¦   ¦   ¦--out_gtf: $( inputs.input_bam.basename.split("_sorted.bam")[0].concat("_exprs.gtf") )   
323  ¦   ¦   ¦--verbose: stringtie_verbose                                                             
324  ¦   ¦   ¦--min_isoform_abundance: stringtie_min_isoform_abundance                                 
325  ¦   ¦   ¦--junction_coverage: stringtie_junction_coverage                                         
326  ¦   ¦   ¦--min_read_coverage: stringtie_min_read_coverage                                         
327  ¦   ¦   ¦--conservative_mode: stringtie_conservative_mode                                         
328  ¦   ¦   ¦--outputdir: $( inputs.input_bam.basename.split("_sorted.bam")[0] )                      
329  ¦   ¦   ¦--expression_estimation_mode: stringtie_expression_estimation_mode                       
330  ¦   ¦   °--ballgown_table_files: stringtie_ballgown_table_files                                   
331  ¦   °--out:                                                                                       
332  ¦       ¦--output_gtf                                                                             
333  ¦       °--outdir                                                                                 
334  ¦--ballgown_de                                                                                    
335  ¦   ¦--run:                                                                                       
336  ¦   ¦   °--../wrappers/ballgown.cwl                                                               
337  ¦   ¦--in:                                                                                        
338  ¦   ¦   ¦--stringtie_dirs: stringtie_expression/outdir                                            
339  ¦   ¦   ¦--phenotype_file: bg_phenotype_file                                                      
340  ¦   ¦   ¦--phenotype: bg_phenotype                                                                
341  ¦   ¦   ¦--samples: bg_samples                                                                    
342  ¦   ¦   ¦--timecourse: bg_timecourse                                                              
343  ¦   ¦   ¦--feature: bg_feature                                                                    
344  ¦   ¦   ¦--measure: bg_measure                                                                    
345  ¦   ¦   ¦--confounders: bg_confounders                                                            
346  ¦   ¦   ¦--custom_model: bg_custom_model                                                          
347  ¦   ¦   ¦--mod: bg_mod                                                                            
348  ¦   ¦   °--mod0: bg_mod0                                                                          
349  ¦   °--out:                                                                                       
350  ¦       ¦--ballgown_de_results                                                                    
351  ¦       ¦--ballgown_object                                                                        
352  ¦       °--ballgown_de_custom_model                                                               
353  ¦--featureCounts                                                                                  
354  ¦   ¦--run:                                                                                       
355  ¦   ¦   °--../wrappers/featureCounts.cwl                                                          
356  ¦   ¦--in:                                                                                        
357  ¦   ¦   ¦--number_of_threads: featureCounts_number_of_threads                                     
358  ¦   ¦   ¦--annotation_file: featureCounts_annotation_file                                         
359  ¦   ¦   ¦--output_file: featureCounts_output_file                                                 
360  ¦   ¦   ¦--inputFiles: samtools_sort/sorted                                                       
361  ¦   ¦   °--read_meta_feature_overlap: featureCounts_read_meta_feature_overlap                     
362  ¦   °--out:                                                                                       
363  ¦       °--output                                                                                 
364  °--DESeq2_analysis                                                                                
365      ¦--run:                                                                                       
366      ¦   °--../wrappers/DESeq2.cwl                                                                 
367      ¦--in:                                                                                        
368      ¦   ¦--count_matrix: featureCounts/output                                                     
369      ¦   ¦--metadata: deseq2_metadata                                                              
370      ¦   ¦--samples: deseq2_samples                                                                
371      ¦   ¦--design: deseq2_design                                                                  
372      ¦   ¦--min_sum_of_reads: deseq2_min_sum_of_reads                                              
373      ¦   ¦--reference_level: deseq2_reference_level                                                
374      ¦   ¦--phenotype: deseq2_phenotype                                                            
375      ¦   ¦--contrast: deseq2_contrast                                                              
376      ¦   ¦--numerator: deseq2_numerator                                                            
377      ¦   ¦--denominator: deseq2_denominator                                                        
378      ¦   ¦--lfcThreshold: deseq2_lfcThreshold                                                      
379      ¦   ¦--pAdjustMethod: deseq2_pAdjustMethod                                                    
380      ¦   ¦--alpha: deseq2_alpha                                                                    
381      ¦   ¦--parallelization: deseq2_parallelization                                                
382      ¦   ¦--cores: deseq2_cores                                                                    
383      ¦   ¦--transformation: deseq2_transformation                                                  
384      ¦   ¦--blind: deseq2_blind                                                                    
385      ¦   ¦--hypothesis: deseq2_hypothesis                                                          
386      ¦   ¦--reduced: deseq2_reduced                                                                
387      ¦   ¦--hidden_batch_effects: deseq2_hidden_batch_effects                                      
388      ¦   ¦--hidden_batch_row_means: deseq2_hidden_batch_row_means                                  
389      ¦   ¦--hidden_batch_method: deseq2_hidden_batch_method                                        
390      ¦   °--variables: deseq2_variables                                                            
391      °--out:                                                                                       
392          ¦--deseq2_de_results                                                                      
393          ¦--deseq2_dds_object                                                                      
394          ¦--deseq2_res_lfcShrink_object                                                            
395          °--deseq2_transformed_object
```
