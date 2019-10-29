# broad-fungalgroup
Broad Fungal Genomics group scripts

This repository hosts scripts used publications to carry out fungal genomic analyses.

## Genotype filtering (filterGatkGenotypes.py)
Usage example:
```
filterGatkGenotypes.py --min_GQ 50 --min_percent_alt_in_AD 0.8 --min_total_DP 10 ${VCF} > ${FILTERED_VCF} 2> variant_qc_genotype_filter.tsv
```
