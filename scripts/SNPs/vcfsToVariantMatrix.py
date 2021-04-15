#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import sys, re
import gzip
import argparse
import vcfTools

parser = argparse.ArgumentParser()
parser.add_argument('infile', help='file of VCF filenames', type=str)
parser.add_argument('--collapse_freq', help='collapse rare variants at a freq <= <freq>, requires annotation', type=float)
parser.add_argument('--max_indel_size', help='ignore insertions > <int> bp', type=int)
parser.add_argument('--ignore_low', help='ignore low impact and intron mutations', action='store_true')
parser.add_argument('--max_amb_samples', help='maximum number of samples with ambiguous calls for a site to be included', type=int)
parser.add_argument('--no_amb_samples', help='do not include any sites with ambiguous calls', action='store_true')
parser.add_argument('--skip_annotations', help='do not print variant effect annotations', action='store_true')
args = parser.parse_args()

infile = args.infile
coll_freq = float(0)
max_indel_size = 0
ignore_low = args.ignore_low
max_amb = 1000000
skip_annotations = args.skip_annotations

if args.collapse_freq:
	coll_freq = args.collapse_freq
if args.max_indel_size:
	max_indel_size = args.max_indel_size
if args.max_amb_samples:
	max_amb = args.max_amb_samples
if args.no_amb_samples:
	max_amb = 0
	
comment_pattern = re.compile(r"^#")

passed_snp_positions = {}
passed_snp_counts = {}
passed_snp_features = {}
amb_snp_positions = {}
amb_snp_counts = {}
annotations = {}
genome_list = []


with open(infile) as fof:
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
						
		if fof_line[-3:] == ".gz":
			vcf_file = gzip.open(fof_line, "rt")
		else:
			vcf_file = open(fof_line)

		for vcf_line in vcf_file:
			if not (re.search(comment_pattern, vcf_line)):
				record = vcfTools.VcfRecord(vcf_line)
				is_passing = record.is_passing(caller)
				for genome in per_file_genome_list:
					chrom = record.get_chrom()
					genotype = record.get_genotype(index=header.get_sample_index(genome),min_gq=0)
					variant_type = record.get_variant_type(caller,genotype)
					alt = record.get_alt(genotype)
					variant_annot = False
					variant_impact = False
					variant_feature = False
					variant_effect = False
					variant_length = 0
					if alt:
						variant_annot = record.get_snpeff_annot(alt)
						variant_impact = record.get_snpeff_impact(variant_annot)
						variant_feature = record.get_snpeff_feature(variant_annot)
						variant_effect = record.get_snpeff_effect(variant_annot)
						variant_length = record.get_variant_length(genotype)	
						
					position = str(record.get_pos())

					## if syn/low and set to ignore syn/low
					## elsif passing snp or passing indel within size range						
					## elsif ambig snp or amb indel within size range
					
					if (ignore_low and (variant_impact == 'LOW' or variant_effect == 'intron_variant')):
						pass
					elif is_passing and (variant_type == 'SNP' or ((variant_type == 'INSERTION' or variant_type == 'DELETION') and (not max_indel_size or variant_length <= max_indel_size))):
						if not genome in passed_snp_positions:
							passed_snp_positions[genome] = {}
						if not chrom in passed_snp_positions[genome]:
							passed_snp_positions[genome][chrom] = {}
						passed_snp_positions[genome][chrom][position] = 1
						
						if not chrom in passed_snp_counts:
							passed_snp_counts[chrom] = {}
							passed_snp_features[chrom] = {}
						if not position in passed_snp_counts[chrom]:
							passed_snp_counts[chrom][position] = 1
							passed_snp_features[chrom][position] = variant_feature
						else:
							passed_snp_counts[chrom][position] += 1
						
						if not chrom in annotations:
							annotations[chrom] = {}
						if not position in annotations[chrom]:
							annotations[chrom][position] = [variant_annot]
						elif not variant_annot in annotations[chrom][position]:
							annotations[chrom][position].append(variant_annot)
							
					elif variant_type == 'uncalled_ambiguous' or (variant_type == 'SNP' or ((variant_type == 'INSERTION' or variant_type == 'DELETION') and (not max_indel_size or variant_length <= max_indel_size))):
						if not genome in amb_snp_positions:
							amb_snp_positions[genome] = {}
						if not chrom in amb_snp_positions[genome]:
							amb_snp_positions[genome][chrom] = {}
						amb_snp_positions[genome][chrom][position] = 1
						
						if not chrom in amb_snp_counts:
							amb_snp_counts[chrom] = {}
						if not position in amb_snp_counts[chrom]:
							amb_snp_counts[chrom][position] = 1
						else:
							amb_snp_counts[chrom][position] += 1
		vcf_file.close()

sys.stderr.write("Collapsing rare variants if enabled...\n")


chroms = sorted(passed_snp_counts.keys())
coll_passed_snp_counts = {}
initial_positions = 0
filtered_positions = 0
final_positions = 0

for chrom in chroms:
	positions = sorted(passed_snp_counts[chrom].keys(), key=int)
	for position in positions:
		initial_positions += 1
		current_amb = 0
		if chrom in amb_snp_counts:
			if position in amb_snp_counts[chrom]:
				current_amb = amb_snp_counts[chrom][position]
		if current_amb <= max_amb:
			filtered_positions += 1
			current_freq = passed_snp_counts[chrom][position]/len(genome_list)
			if current_freq <= coll_freq:
				current_feature = passed_snp_features[chrom][position]
				if current_feature:
					if not chrom in coll_passed_snp_counts:
						coll_passed_snp_counts[chrom] = {}
					if not current_feature in coll_passed_snp_counts[chrom]:
						coll_passed_snp_counts[chrom][current_feature] = 1
						final_positions += 1
					
					current_annots = annotations[chrom][position]
					if not current_feature in annotations[chrom]:
						annotations[chrom][current_feature] = current_annots
					else:
						annotations[chrom][current_feature].extend(current_annots)
				
					for genome in genome_list:
						if genome in passed_snp_positions:
							if chrom in passed_snp_positions[genome]:
								if position in passed_snp_positions[genome][chrom]:
									passed_snp_positions[genome][chrom][current_feature] = 1
								elif genome in amb_snp_positions and chrom in amb_snp_positions[genome] and position in amb_snp_positions[genome][chrom]:
									amb_snp_positions[genome][chrom][current_feature] = 1
				else:
					sys.stderr.write("Cannot collapse " + chrom + " position " + str(position) + ": no annotated feature. Skipping...\n")
			else:
				if not chrom in coll_passed_snp_counts:
					coll_passed_snp_counts[chrom] = {}
				coll_passed_snp_counts[chrom][position] = 1
				final_positions += 1

sys.stderr.write(str(initial_positions) + " initial positions, filtered to " + str(filtered_positions) + " positions, collapsed to " + str(final_positions) + " final positions.\n")

print("#feature", end="")
for chrom in chroms:
	positions = sorted(coll_passed_snp_counts[chrom].keys())
	for position in positions:
		print("".join(["\t",chrom,"_",str(position)]), end="")
print("")

if not skip_annotations:
	print("#annotations", end="")
	for chrom in chroms:
		positions = sorted(coll_passed_snp_counts[chrom].keys())
		for position in positions:
			print("\t", end="")
			try:
				print(",".join(annotations[chrom][position]), end="")
			except:
				print("NONE", end="")
	print("")

print("reference", end="")
for chrom in chroms:
	positions = sorted(coll_passed_snp_counts[chrom].keys())
	for position in positions:
		print("\t0", end="")
print("")

for genome in genome_list:
	print(genome, end="")
	for chrom in chroms:
		positions = sorted(coll_passed_snp_counts[chrom].keys())
		for position in positions:
			if genome in passed_snp_positions:
				if chrom in passed_snp_positions[genome]:
					if position in passed_snp_positions[genome][chrom]:
						print("\t1", end="")
					elif genome in amb_snp_positions and chrom in amb_snp_positions[genome] and position in amb_snp_positions[genome][chrom]:
						print("\t-", end="")
					else:
						print("\t0", end="")
	print("")
	


		
