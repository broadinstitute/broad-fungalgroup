#!/usr/bin/env python

import sys, re
import argparse
import vcfTools, gffReader
from Bio import SeqIO
from Bio.Seq import Seq



parser = argparse.ArgumentParser()
parser.add_argument('vcf_list_file', help='file of VCF filenames', type=str)
parser.add_argument('fasta_file', help='reference fasta', type=str)
parser.add_argument('gff_file', help='reference gff', type=str)
parser.add_argument('gene', help='gene', type=str)
parser.add_argument('--include_ref', help='include reference in the matrix', action='store_true')
args = parser.parse_args()

infile = args.vcf_list_file
fasta_file = args.fasta_file
gff_file = args.gff_file
gene = args.gene
include_ref =args.include_ref

	
ref_db = {}
	
fasta = open(fasta_file, "rU")
for record in SeqIO.parse(fasta, "fasta"):
    ref_db[record.id] = record.seq
fasta.close()

#print(sequences)

gff = gffReader.GffRecord(gff_file)

positions = gff.get_gene_positions(gene)
gene_chrom = gff.get_gene_seq_id(gene)
gene_strand = gff.get_gene_strand(gene)

#print(positions)


comment_pattern = re.compile(r"^#")

ref_bases = {}
alt_bases = {}
passed_snp_positions = {}
genome_list = []
amb_pos_counts = {}

with open(infile, 'r') as fof:
	for fof_line in fof:
		sys.stderr.write("Searching " + fof_line)
		fof_line = fof_line.rstrip()
		header = vcfTools.VcfHeader(fof_line)
		caller = header.get_caller()
		samples = header.get_samples()
		per_file_genome_list = []
		for sample in samples:
			genome_list.append(sample)
			per_file_genome_list.append(sample)
						
		with open(fof_line, 'r') as vcf_file:
			for vcf_line in vcf_file:
				if not (re.search(comment_pattern, vcf_line)):
					record = vcfTools.VcfRecord(vcf_line)
 					chrom = record.get_chrom()
 					break_flag = False
					if chrom == gene_chrom:
	 					pos = int(record.get_pos())
	 					#print(pos)
	 					if pos in positions:
	 						#print("yes!")
	 						break_flag = True	
							pass_or_fail = record.is_passing(caller)
							for genome in per_file_genome_list:
								genotype = record.get_genotype(index=header.get_sample_index(genome),min_gq=0)
								variant_type = record.get_variant_type(caller,genotype)																		
								if pass_or_fail and not variant_type:					
									pass
								elif pass_or_fail and variant_type == 'SNP':
									chrom = record.get_chrom()
									pos = int(record.get_pos())
					
									if not genome in alt_bases:
										alt_bases[genome] = {}
									if not chrom in alt_bases[genome]:
										alt_bases[genome][chrom] = {}
									alt_bases[genome][chrom][pos] = record.get_alt(genotype)
									### print(alt_bases[genome][chrom][pos]) ###
									### print(alt_bases) ###
					
									if not chrom in ref_bases:
										ref_bases[chrom] = {}						
									ref_bases[chrom][pos] = record.get_ref()
					
									if not chrom in passed_snp_positions:
										passed_snp_positions[chrom] = {}
									passed_snp_positions[chrom][pos] = True
								else:
									chrom = record.get_chrom()
									pos = int(record.get_pos())
									if not genome in alt_bases:
										alt_bases[genome] = {} 
									if not chrom in alt_bases[genome]: 
										alt_bases[genome][chrom] = {}							
									alt_bases[genome][chrom][pos] = 'N'
						elif break_flag == True:
							break


#### dont forget to account for RC  ######

if include_ref:
	genome_list.append('reference')

for ind_genome in genome_list:
	print(">" + ind_genome)
	sequence = ''
	for position in positions:
		try:
			if position in alt_bases[ind_genome][gene_chrom]:
				sequence += alt_bases[ind_genome][gene_chrom][position]
			else:
				sequence += ref_bases[gene_chrom][position]
		except:
			try:
				sequence += ref_bases[gene_chrom][position]
			except:
				sequence += ref_db[gene_chrom][(position-1)]
				#print(ref_db[gene_chrom].seq)
	if gene_strand == '-':
		sequence = Seq(sequence)
		sequence = sequence.reverse_complement()
	for i in range(0,len(sequence),60):
		print(sequence[i:i+60])
