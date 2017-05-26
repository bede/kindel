---
title: 'Kindel: indel-aware consensus for nucleotide sequence alignments'
tags:
  - bioinformatics
  - sequence analysis
  - genome assembly 
authors:
 - name: Bede Constantinides
   orcid: 0000-0002-3480-3819
   affiliation: 1
 - name: David L. Robertson
   orcid: 0000-0001-6338-0221
   affiliation: 2
affiliations:
 - name: Evolution and Genomic Sciences, University of Manchester, Manchester, UK
   index: 1
 - name: MRC-University of Glasgow Centre for Virus Research, Glasgow, UK
   index: 2
date: 26 May 2017
---



# Summary

Kindel is a collection of command line tools built around a Python package for inferring consensus sequence from an alignment of nucleotide sequences in Sequence Alignment/Map (SAM) format (Li *et al.* 2009). Sequence substitutions, insertions and deletions are represented in the generated consensus where the majority of  sequence information—even that contained within unaligned sequence regions—supports their presence. In this way, Kindel  generates a data-specific reference sequence with maximum similarity to mapped sequences for use in phylogenetic and other analyses. 

Kindel may additionally be used both to scrutinise and visualise low frequency variation present in diverse populations through generation of positionwise feature tables, facilitating, for example, the comparison of intrapatient viral populations sampled from multiple individuals and/or timepoints.



# References

**Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R** (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics,* 25 (16): 2078-2079, doi:10.1093/bioinformatics/btp352