# broad-fungalgroup
Broad Fungal Genomics group scripts

This repository hosts scripts used publications to carry out fungal genomic analyses.

## Genotype filtering (filterGatkGenotypes.py)
Usage example:
```
filterGatkGenotypes.py --min_GQ 50 --min_percent_alt_in_AD 0.8 --min_total_DP 10 ${VCF} > ${FILTERED_VCF} 2> variant_qc_genotype_filter.tsv
```

## create whole genome fasta-alignment from VCF (SNPs only) (vcfSnpsToFasta.py)
Usage example:
```
list_vcf = file containing the name of the input VCF
python vcfSnpsToFasta.py --max_amb_samples 40 list_vcf > filtered.SNPs.fa
```

## create gene fasta-alignment from VCF (SNPs only) (vcfSnpsToGeneFasta.py)
Usage example:
```
python vcfSnpsToGeneFasta.py list_of_vcf_files reference.fasta reference.gff gene_id > output.fasta
```
