# AsGARD
## *An. stephensi* Genomics and Resistance Diagnostics

This repo will contains analysis code supporting the CEASE (https://www.ceasemalaria.org) genomics work package.

Tristan Dennis, 2024

## The data.

We have generated a mature catalogue of single nucleotide polymorphism (SNP), phased haplotype, and copy number variant (CNV) calls, for 475 _An. stephensi_ samples. In addition, we have generated an accessibility mask and draft recombination map. This genomic toolkit will be made freely accessible to to download and analyse. We are still finalising the particulars of storage and availability, but more information will be available soon. The raw reads have been uploaded to ENA but are embargoed until we submit our preprint.

## Variant discovery and format.

The raw read data were used to generate SNP calls using NVIDIA's Clara Parabricks [Germline](https://docs.nvidia.com/clara/parabricks/3.7.0/documentation/tooldocs/man_germline.html) workflow. This sped up our per-sample read-variant calling to around 11 minutes per sample, and helps get around issues that other groups have had previously using HaplotypeCaller in highly diverse species like mosquitoes. The data were then filtered, mostly according to the regimen used by [MalariaGen](https://malariagen.github.io/vector-data/ag3/methods.html), and converted to Zarr format. 

[Zarr](https://zarr.readthedocs.io/en/stable/) is a format for compressed, chunked, n-dimensional arrays, that facilitates analysis of genomic data quickly and efficiently (an example is that you don't have to load the whole dataset into memory to perform an analysis). The scripts for doing so are in AsGARD/scripts. 

## CNV calling

The CNV calling was performed using [this Hidden Markov Model](https://pubmed.ncbi.nlm.nih.gov/27531718/), according to the method published in [Eric's 2019 paper](https://genome.cshlp.org/content/26/9/1288.long). The original code is located in [this repo](https://github.com/EricRLucas/CNV_pipeline). We have vague plans to tidy it up into a WDL/Nextflow at some point. The modified scripts are in this repository.

## Phasing

The data were phased using BEAGLE. See the manuscript or code for exact parameters.

## Analysis

Most of the analysis in this repository uses the incomparable [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/). Even though AsGARD is not part of MalariaGen, I work with the MalariaGen team and, as exemplars of vector data analysis, a lot of the analysis code will be familiar to malariagen-data developers and users ;). 

The code required to reproduce the analyses in [THE PAPER] is located in AsGARD/analysis. The directory is organised into various subfolders, and is a mixture of Jupyter Notebooks, R, shell, and python scripts. For maps and stats, I've mainly used R, for scripting, Python and Shell and for population genetics, Python.

## Data availability

When the manuscript is published, SNP and CNV data will be available in MalariaGEN Data.



