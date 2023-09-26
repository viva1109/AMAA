# AMAA: All-inclusive Microbiome Association Analysis

## Overview
Preprocessing steps such as OTU clustering, choice of OTU filtering threshold, rarefying can be a source of heterogeneity. AMAA is a package that envelope a pipeline building microbial count tables based on different databases and clustering method and the methods for metagenome-wide association analysis. It provides the convenient use of various methods for microbiome association analysis under a unified preprocessing and comparison the results based on different databases or clustering method.
### Mainteinor
Kangjin Kim, Ph.D Candidate <rekki@channing.harvard.edu>
Sungho Won, Ph.D <won1@snu.ac.kr>

## Installation
### Get tools and example datsets from github
    git clone https://github.com/viva1109/amaa
### Go to cloned folder (default: ./amaa)
    cd amaa
### Download and install Docker
Docker is available at https://docker.com
### Get and run AMAA docker image
    docker run --name amaa -v $(pwd):/home/ -it viva1109/amaa:1.02
## Example 1: Sequencing data processing
### Analysis manual
    bash 16S_pipeline/pipeline_asv.sh -h
### Example run
    bash 16S_pipeline/pipeline_asv.sh -d /home/examples/AMAAtoy/fastq -c /bash home/16S_pipeline/pip.configure -z silva -p open -w T -s 1 -e 11 -o /home/examples/AMAAtoy/table -t 1 -f CCTACGGGNGGCWGCAG -r GACTACHVGGGTATCTAATCC
## Example 2: Cross-sectional metagenome-wide association analysis
### Analysis manual
    bash 16S_pipeline/Stat_Analysis_pipeline.sh -h
### Example run
    bash 16S_pipeline/Stat_Analysis_pipeline.sh -d /home/examples/CrossSectional/data/ -D /home/examples/CrossSectional/result -t 2 -z silva -s 1 -e 3 -M TMAT15_oMirkat_ANCOM -p 1000 -c 3 -P pheno_togo.csv -C 0 -n Control_Case
## Example 3: Longitudinal metagenome-wide association analysis
### Analysis manual
    bash 16S_pipeline/Stat_Analysis_pipeline_longitudinal.sh -h
### Example run
    bash 16S_pipeline/Stat_Analysis_pipeline_longitudinal.sh -d /home/examples/Longitudinal/data/ -D /home/examples/Longitudinal/result -t 2 -z ez -M LMM_NBMM_ZINBMM_mTMAT -P pheno_togo.csv -C 0 -n Control_Case -V BinTrait -v Cov1,Time -Z NULL -T Time

