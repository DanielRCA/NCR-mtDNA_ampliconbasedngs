# NCR-mtDNA_ampliconbasedngs

## Introduction

DanielRCA/NCR-mtDNA_ampliconbasedngs is an easy-to-use bioinformatic tool to analyse the [Non-Coding Region](https://en.wikipedia.org/wiki/MtDNA_control_region) of ancient human mtDNA obtained by an amplicon-based Next Generation Sequencing method ([PowerSeq<sup>(TM)</sup> CRM Nested System kit](https://www.promega.es/products/forensic-dna-analysis-mps/target-amplification-and-library-prep/powerseq-crm-nested-system-custom/?catNum=AX5810), Promega Corporation).

## Pipeline Summary

The pipeline performs the following:

* Reference genome indices creation for mapping (`bwa` and `samtools`)
* Sequencing quality control (`FastQC`)
* Sequencing adapter removal, duplicates removal, primers removal, paired-end data merging, Illumina two-coloured sequencer poly-G tail removal, post-adapter-removal trimmnig of FASTQ files prior mapping (`fastp`)
* Read mapping to reference (`bwa aln` and `bwa sampe`)
* Post-mapping processing, statistics and conversion to bam (`samtools`)
* Damaged reads extraction (`PMDtools`)
* Post-mapping statistics and BAM quality control (`Qualimap`)
* Creation of VCF genotyping files (`freebayes` and `vcflib`)
* Haplogroup determination (`HaploGrep2`)


## Pre-requesites

A group of tools must be pre-installed. For each tool, version used by our group is shown in brackets. We provide a file with a Conda's environment with all the tools already installed in order to make easier the installation (see [**Installation and usage**](##installation-and-usage)).

- FastQC (v0.11.9)
- fastp (v0.23.2)
- BWA (v0.7.17)
- SAMtools (v1.16.1)
- QualiMap (v2.2.2a)
- PMDtools (v0.60)
- freebayes (v1.3.6)
- vcflib (v1.0.3)
- Haplogrep (v2.4.0)

> [!WARNING]
> The use of other versions, especially concerning freebayes and the VCF generation, could modify the results.

##Installation and usage

You may need a [Linux Operating System](https://en.wikipedia.org/wiki/Linux). You should also have [Conda](https://docs.conda.io/projects/conda/en/stable/) installed, which can be done by following [Anaconda's User Guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). If you are not accustomed to using Linux, we strongly recommend following [this tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-the-anaconda-python-distribution-on-ubuntu-22-04) created in 2022 by Lisa Tagliaferri and Tony Tran.

Once Anaconda has been installed, open the [Command-line interface](https://en.wikipedia.org/wiki/Command-line_interface#Command_prompt) by typing `command line`, `cmd`, or `prompt` into the search bar of your computer.

...

### Inputs

...

### Outputs

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
