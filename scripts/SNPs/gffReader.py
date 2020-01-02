# python vcf tools

from __future__ import print_function
from __future__ import division

import re

class GffRecord:
	"""
	
	"""
	def __init__(self, gff_file):
		self.strand = {}
		self.gene_to_transcript = {}
		self.gene_to_strand = {}
		self.gene_to_sequence = {}
		self.gene_to_name = {}
		self.transcript_to_cds_coords = {}
		self.ends_to_genes = {}   ### right now for single scaffold only
		self.sorted_ends = []
		self.last_end = 0  ### right now for single scaffold only
		
		with open(gff_file, 'r') as gff:
			for line in gff:
				fields = line.rstrip().split('\t')
				try:
					if fields[2] == 'mRNA':
						id = re.search('ID=([^;]+)',fields[8])
						parent = re.search('Parent=([^;]+)',fields[8])
						self.gene_to_transcript[parent.group(1)] = id.group(1)
						self.gene_to_strand[parent.group(1)] = fields[6]
						self.gene_to_sequence[parent.group(1)] = fields[0]
				except:
					pass
		
				try:
					if fields[2] == 'CDS':
						parent = re.search('Parent=([^;]+)',fields[8])
						if not parent.group(1) in self.transcript_to_cds_coords:
							self.transcript_to_cds_coords[parent.group(1)] = []
						self.transcript_to_cds_coords[parent.group(1)].append((fields[3],fields[4]))	
				except:
					pass
					
				try:
					if fields[2] == 'gene':   ### right now for single scaffold only
						id = re.search('ID=([^;]+)',fields[8])
						name = re.search('Name=([^;]+)',fields[8])
						self.gene_to_name[id.group(1)] = name.group(1)
						self.ends_to_genes[int(fields[3])] = id.group(1)
						self.ends_to_genes[int(fields[4])] = id.group(1)
						if int(fields[4]) > self.last_end:
							self.last_end = int(fields[4])
				except:
					pass
		
		self.sorted_ends = sorted(self.ends_to_genes.keys())
		##print(self.sorted_ends)
		##print(self.gene_to_name)
						
					
	def get_gene_positions(self, gene):
		positions = []
		transcript = self.gene_to_transcript[gene]
		cds_tuples = self.transcript_to_cds_coords[transcript]
		for coord in cds_tuples:
			for i in range(int(coord[0]),int(coord[1])+1):
				positions.append(i)
		positions.sort()
		return positions
		
	def get_gene_seq_id(self, gene):
		seq_id = self.gene_to_sequence[gene]
		return seq_id
		
	def get_gene_strand(self, gene):
		strand = self.gene_to_strand[gene]
		return strand
		
	def get_gene_name(self,gene):
		try:
			name = self.gene_to_name[gene]
			return name
		except:
			return("None")
		
	def get_feature_from_position(self, position):  ### right now for single scaffold only
		left_end = 0
		right_end = 1000000000
		left_feature = 0
		right_feature = 0
		self.gene_to_strand[0] = '+'
		position = int(position)
		for sorted_end in self.sorted_ends:
			if position > sorted_end and sorted_end > left_end:
				left_end = sorted_end
				left_feature = self.ends_to_genes[sorted_end]
			if position < sorted_end and sorted_end < right_end:
				right_end = sorted_end
				right_feature = self.ends_to_genes[sorted_end]
		if left_feature == right_feature:
			return left_feature
		else:
			combined_feature = "_".join([str(left_feature),self.gene_to_strand[left_feature],str(right_feature),self.gene_to_strand[right_feature]])
			return combined_feature
				
	def get_feature_name_from_position(self, position):  ### right now for single scaffold only
		left_end = 0
		right_end = 1000000000
		left_feature = 0
		right_feature = 0
		self.gene_to_strand[0] = '+'
		position = int(position)
		for sorted_end in self.sorted_ends:
			if position > sorted_end and sorted_end > left_end:
				left_end = sorted_end
				left_feature = self.ends_to_genes[sorted_end]
			if position < sorted_end and sorted_end < right_end:
				right_end = sorted_end
				right_feature = self.ends_to_genes[sorted_end]
		if left_feature == right_feature:
			name = self.get_gene_name(left_feature)
			return name
		else:
			lname = self.get_gene_name(left_feature)
			rname = self.get_gene_name(right_feature)
			combined_name = ",".join([str(lname),str(rname)])
			return combined_name	
	