# AsGARD
## An. stephensi Genomics and Resistance Diagnostics

This repo will contain all the code for the CEASE (https://www.ceasemalaria.org) genomics work package. At the moment, it contains a few helper scripts and analysis notebooks (notebooks), as I tidy everything off of LSTM's servers to prepare these data for publication.

Tristan Dennis, 2024

## The data.

We have generated a mature catalogue of single nucleotide polymorphism (SNP), phased haplotype, and copy number variant (CNV) calls, for 475 _An. stephensi_ samples. In addition, we have generated an accessibility mask and draft recombination map. This genomic toolkit will be made freely accessible to to download and analyse. We are still finalising the particulars of storage and availability, but more information will be available soon. The raw reads have been uploaded to ENA but are embargoed until we submit our preprint.

## Variant discovery and format.

The raw read data were used to generate SNP calls using NVIDIA's Clara Parabricks [Germline](https://docs.nvidia.com/clara/parabricks/3.7.0/documentation/tooldocs/man_germline.html) workflow. This sped up our per-sample read-variant calling to around 11 minutes per sample, and helps get around issues that other groups have had previously using HaplotypeCaller in highly diverse species like mosquitoes. The data were then filtered, mostly according to the regimen used by [MalariaGen](https://malariagen.github.io/vector-data/ag3/methods.html), and converted to Zarr format. 

[Zarr](https://zarr.readthedocs.io/en/stable/) is a format for compressed, chunked, n-dimensional arrays, that facilitates analysis of genomic data quickly and efficiently (an example is that you don't have to load the whole dataset into memory to perform an analysis). The scripts for doing so are in AsGARD/scripts. 

## Analysis

Most of the analysis in this repository uses the incomparable [scikit-allel](https://scikit-allel.readthedocs.io/en/stable/). Even though AsGARD is not part of MalariaGen, I work with the MalariaGen team and, as exemplars of vector data analysis, a lot of the analysis code will be familiar to malariagen-data developers and users ;). I am steadily tidying up my notebooks and making them available in AsGARD/notebooks.
