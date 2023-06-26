# CWL pipelines for the analysis of Next-Generation Sequencing data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description
Here, we established automated software pipelines for the asychronous analysis of Next-Generation Sequencing (NGS), single- and/or paired-end data from: 

- RNA-Seq
- ChIP-Seq
- Germline variant calling 

experiments by incorporating various software tools into functional workflows using the [Common Workflow Language](https://www.commonwl.org/) (CWL v1.0), an open standard for describing computational workflows for data-intensive science fields such as Bioinformatics ([Amstutz et al., 2016](https://doi.org/10.6084/m9.figshare.3115156.v2)). The pipelines are implemented as [CWL Workflows](https://www.commonwl.org/v1.2/Workflow.html) and include a number of wrappers, in the form of [CWL CommandLineTools](https://www.commonwl.org/v1.0/CommandLineTool.html), for various computational biology tools (Trimmomatic, HISAT2, BWA-MEM, samtools, etc.). To resolve any issues regarding software dependencies and compatibility, the ability of CWL to support operations through Docker containers was utilized. Docker is a containerization platform that allows for packaging an application along with its dependencies and running it in a self-contained unit called a container ([Merkel, 2014](https://www.linuxjournal.com/content/docker-lightweight-linux-containers-consistent-development-and-deployment)).

## Setup

### Requirements

The workflow requires a [CWL](https://www.commonwl.org/) runner and [Docker](https://docs.docker.com/) to operate. It is recommended to follow the detailed installation instructions available in [CWL tool description reference implementation](https://github.com/common-workflow-language/cwltool). For those with no experience in CWL, the [CWL user guide](https://www.commonwl.org/user_guide/) is recommended as a starting point.

### Installation

Clone the git repository:

```bash
# With https (default)
git clone https://github.com/BiodataAnalysisGroup/cwl_pipelines.git
```
```bash
# With ssh
git clone git@github.com:BiodataAnalysisGroup/cwl_pipelines.git
```

### Usage

Configuration of inputs and parameters is performed using a YAML file, with an example template, that is thoroughly commented to simplify the process, being available for every workflow in `<workflow-of-interest>/yaml_files`.

```bash
# use CWL runner (e.g., cwltool)
cwltool <workflow-of-interest>/workflows/<workflow-of-interest>.cwl <workflow-of-interest>/yaml_files/<workflow-of-interest>.yml
```
An easy way to keep a log for a running workflow includes a few extra shell commands:

```bash
#!/bin/bash
# save stdout and stderr in log file and also show on screen
exec > >(tee -a <log-file>.out) 2>&1
# mark starting date and time in log file:
echo $(date +%Y_%m_%d-%T)
# set custom directory and prefix for temporary (tmp) files
prefix=/<tmp-directory>/<tmp-prefix>
# e.g., prefix=/data/TMP_directory/TMP_ 
# execute workflow
cwltool --parallel --rm-container --rm-tmpdir --tmpdir-prefix $prefix <workflow-of-interest>/workflows/<workflow-of-interest>.cwl <workflow-of-interest>/yaml_files/<workflow-of-interest>.yml
# mark ending date and time in log file:
echo $(date +%Y_%m_%d-%T)
```
### Input files
The worklflow accepts FASTQ files, `.gz` compressed or not. The user is required to specify the directory where the files reside. The workflow will automatically detected files with `.fq(.gz)`/`.fastq(.gz)` suffix as FASTQ format inputs, respectively. Furthermore, a file name convention is used in order to automatically identify input files, separate (i) single- and paired-end files as well as (ii) forward and reverse reads, and assign sample identifiers for downstream analysis. The file name convention is shown below:

<p align="center">
<img src="https://github.com/BiodataAnalysisGroup/kmerCountClassifier/blob/main/file_name_convention.png" alt="file name convention" width="500">
</p>

where ``sample_name`` corresponds to the unique file name of each sample, ``file_split`` defines the separator string pattern between the sampleName and the rest of the file name, ``fwd_pattern`` and ``rev_pattern`` the string patterns that define single-end or paired-end forward (contains ``fwd_pattern``) and paired-end reverse files (contains ``rev_pattern``). For example, the single-end FASTQ file **SAMN18116076_1.fastq.gz** contains SAMN18116076 as unique ``sample_name``, `_` as file_split and `_1` as fwd corresponding to a single-end file. It is the same for the files **SAMEA3751395_1.fastq.gz** and **SAMEA3751395_2.fastq.gz**, except that here two files have the same ``sample_name``, with one containing the ``fwd_pattern`` and the other the ``rev_pattern``, thus identified as paired-end files of the same sample. 

These parameters can be specified within the YAML template, with ``file_split``, ``fwd_pattern`` and ``rev_pattern`` corresponding to the ``input_file_split``, ``input_file_split_fwd_single`` and ``input_file_split_rev`` parameters, respectively.

## References

- Chapman, B. et al. (2016) Common Workflow Language, v1.0. Edited by P. Amstutz, M.R. Crusoe, and N. Tijanić. United States: figshare. Available at: https://doi.org/10.6084/m9.figshare.3115156.v2.
- Merkel, D. (2014) ‘Docker: lightweight linux containers for consistent development and deployment’, Linux journal, 2014(239), p. 2.
