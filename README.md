# NCR-mtDNA_ampliconbasedngs

## Introduction

DanielRCA/NCR-mtDNA_ampliconbasedngs is an easy-to-use bioinformatic tool to analyse the [Non-Coding Region](https://en.wikipedia.org/wiki/MtDNA_control_region) of ancient human mtDNA obtained by an [amplicon](https://en.wikipedia.org/wiki/Amplicon)-based [Next Generation Sequencing](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3841808/) method ([PowerSeq<sup>(TM)</sup> CRM Nested System kit](https://www.promega.es/products/forensic-dna-analysis-mps/target-amplification-and-library-prep/powerseq-crm-nested-system-custom/?catNum=AX5810), Promega Corporation).

## Pipeline Summary

The pipeline performs the following:

1. Reference genome indices creation for mapping (`bwa` and `samtools`)
2. Sequencing quality control (`FastQC`)
3. Sequencing adapter removal, duplicates removal, primers removal, paired-end data merging, Illumina two-coloured sequencer poly-G tail removal, post-adapter-removal trimmnig of FASTQ files prior mapping (`fastp`)
4. Read mapping to reference (`bwa aln` and `bwa sampe`)
5. Post-mapping processing, statistics and conversion to bam (`samtools`)
6. Damaged reads extraction (`PMDtools`)
7. Post-mapping statistics and BAM quality control (`Qualimap`)
8. Creation of VCF genotyping files (`freebayes` and `vcflib`)
9. Haplogroup determination (`HaploGrep2`)


## Pre-requesites

A group of tools must be pre-installed. For each tool, version used by our group is shown in brackets. We provide a Conda enviroment file with all the versions indicated to make installation easier installation (see [Installation and usage](#installation-and-usage)).

- FastQC (v0.11.9)
- fastp (v0.23.2)
- BWA (v0.7.17)
- SAMtools (v1.16.1)
- QualiMap (v2.2.2a)
- PMDtools (v0.60)
- freebayes (v1.3.6)
- vcflib (v1.0.3)
- Haplogrep (v2.4.0)

> [!IMPORTANT]
> The use of other versions, especially  in the context of freebayes and VCF generation, could alter results.

## Installation and usage

### Installation

You will need a [Linux Operating System](https://en.wikipedia.org/wiki/Linux). You must also have [Conda](https://docs.conda.io/projects/conda/en/stable/) installed, which can be done by following [Anaconda's User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). If you are not accustomed to using Linux, we strongly recommend following [this tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-the-anaconda-python-distribution-on-ubuntu-22-04) created in 2022 by Lisa Tagliaferri and Tony Tran.

Once Anaconda has been installed, download [Conda_Environment_File_NCR-mtDNA_ngsamplicon.txt](https://github.com/DanielRCA/NCR-mtDNA_ampliconbasedngs/blob/main/Conda_Environment_File_NCR-mtDNA_ngsamplicon.txt) file and save it in a directory. This is the [Conda environement](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html) file. Next, open the [Command-line interface](https://en.wikipedia.org/wiki/Command-line_interface#Command_prompt) by typing `Ctrl+Alt+T` in the directory where you saved your Conda environment file. Alternatively, you can type `command line`, `cmd`, or `prompt` into the search bar of your computer and navigate to the directory where you save the Conda environment file using `cd` command. If you face problems opening the command-line, here is a [good tutorial](https://www.groovypost.com/howto/cant-open-terminal-in-ubuntu-fixes/).

Install the environment with the command:
  ```bash
  conda create --name env_ncrngsamplicons --file Conda_Environment_File_NCR-mtDNA_ngsamplicon.txt
  ````
You can open the environment you have just created with the command:
  ```bash
  conda activate env_ncrngsamplicons
  ````
And close it typing:
  ```bash
  conda deactivate
  ````
You can chanage the environment's name: `env_ncrngsamplicons`.

### Usage

You have to be in the directory where you have all your input files and type on the prompt:
  ```bash
  bash NCR_ampliconbasedngs.sh
  ```
We strongly recommend to save the information in a log file, so you can find where the possible mistake is (time will give you how much time the script was running and can ve avoid). Type as follow:
```bash
  time bash NCR_ampliconbasedngs.sh 2>&1 | tee NCR_ampliconbasedngs.log
  ```

#### Inputs


NCR_ampliconbasedngs.sh, rCRS_NCR_lineal.fasta and range.py must be downloaded and saved in a directory

...

> [!WARNING]
> If any sample has the word '**tem**' in its name, some problems could appear since this combination of letters is used to name temporary files that would be removed at the end of the pipeline. If it is in capital letters '**TEM**', no problem should occur.

#### Outputs

Once the indices for the reference genome are created, they will be saved so that they do not have to be created each time.

Several outputs are generated:

 - A table containing the following information: ID, number of initial reads, duplication rate, number of useful reads, percentage of useful reads, mean depth coverage, mean mapping quality, number of useful damaged reads, percentage of damaged reads, number of mixed bases, haplogroup, haplogroup quality, range of positions and haplotype, all for depth coverage ≥ 5 and depth coverage ≥ 10.
 - An HSD file to upload to Haplogrep3, wich provides more specific information not available through local analysis. This file contains all the information for both all reads and only reads with post-mortem molecular damage extracted by PMDtools for depth coverage ≥ 5 and ≥ 10.
 - Two TXT files, one for each depth coverage analysed, containing information on the mixed bases (sample, position, alternative base, depth coverage of the position, number of appearances of the alternative base, and percentage of the alternative base). A base only will appear if the percentage of mixture is between 30 and 70%.
  
![pipeline drawio](https://github.com/DanielRCA/NCR-mtDNA_ampliconbasedngs/assets/97441691/a9366d58-d987-4771-a6ed-ea2622bc18cb)

## Citations

This pipeline is under revision. Please, cite this repository.

References of the tools used:
- **FastQC** S. Andrews (2010); FASTQC. A quality control tool for high throughput sequence data, https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **fastp** S. Chen, Y. Zhou, Y. Chen, J. Gu (2018), fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, Volume 34, Issue 17, pages i884–i890. https://doi.org/10.1093/bioinformatics/bty560
- **BWA** H. Li, R. Durbin, (2009); Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, Volume 5, pages 1754–1760. https://doi.org/10.1093/bioinformatics/btp324
- **SAMtools** P. Danecek, JK. Bonfield, J. Liddle, J. Marshall, V. Ohan, MO. Pollard, A. Whitwham, T. Keane, SA. McCarthy, RM. Davies, H. Li (2021); Twelve years of SAMtools and BCFtools. GigaScience, Volume 10, Issue 2. https://doi.org/10.1093/gigascience/giab008
- **QualiMap** K. Okonechnikov, A. Conesa, F. García-Alcalde (2016); Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data.  Bioinformatics, Volume 32, Issue 2, pages 292–294. https://doi.org/10.1093/bioinformatics/btv566
- **PMDtools** P. Skoglund, BH. Northoff, MV. Shunkov, A. Derevianko, S. Pääbo, J. Krause, M. Jakobsson (2014); Separating ancient DNA from modern contamination in a Siberian Neandertal. Proceedings of the National Academy of Sciences USA. https://doi.org/10.1073/pnas.1318934111
- **freebayes** E. Garrison, G. Marth (2012); Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN]
- **vcflib** E. Garrison, ZN. Kronenberg, ET. Dawson, BS Pedersen, P Prins (2022); A spectrum of free software tools for processing the VCF variant call format: vcflib, bio-vcf, cyvcf2, hts-nim and slivar. PLoS Comput Biol Volume 18, Issue 5. https://doi.org/10.1371/journal.pcbi.1009123
- **HaploGrep2** H. Weissensteiner, D. Pacher, A. Kloss-Brandstätter, L. Forer, G. Specht, HJ. Bandelt, F. Kronenberg, A. Salas, S. Schönherr (2016); HaploGrep 2: mitochondrial haplogroup classification in the era of high-throughput sequencing. Nucleic Acids Res. Volume 44, pages 58–63. https://doi.org/10.1093/nar/gkw233
