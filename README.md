# NCR-mtDNA_ampliconbasedngs

## Introduction

DanielRCA/NCR-mtDNA_ampliconbasedngs is an easy-to-use bioinformatic tool to analyse the Non-Coding Region of ancient human mtDNA obtained by an amplicon-based Next Generation Sequencing method (PowerSeqTM CRM Nested System kit, Promega Corporation).

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

A group of tools must be pre-installed. For each tool, version used by our group is shown in brackets. The use of other versions could modify the results. We provide a file with a Conda's environment with all the tools already installed in order to make easier the installation (see **Installation and usage**).

- FastQC (v0.11.9)
- fastp (v0.23.2)
- BWA (v0.7.17)
- SAMtools (v1.16.1)
- QualiMap (v2.2.2a)
- PMDtools (v0.60)
- freebayes (v1.3.6)
- vcflib (v1.0.3)
- Haplogrep (v2.4.0)


## Installation and usage

You may need a Unix OS and you should have [Conda](https://docs.conda.io/projects/conda/en/stable/) by installing **Anaconda** follwong the [User guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). If you are not use to Linux, we strongly recommend follow [this tutorial](https://www.digitalocean.com/community/tutorials/how-to-install-the-anaconda-python-distribution-on-ubuntu-22-04) created in 2022 by Lisa Tagliaferri and Tony Tran.

### Inputs

### Outputs

## Citations


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
