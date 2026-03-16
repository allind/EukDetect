#! /usr/bin/env python
#usage: script.py [bamfile] [buscos]
from pysam import AlignmentFile
from Bio import SeqIO
import sys

def main(argv):

	#parse reference file..?
	#parse samtools header
	with AlignmentFile(argv[1]) as bam:
		# Look up time O(1) instead of O(n) with list
		observed = dict.fromkeys(read.reference_name for read in bam.fetch())

		ref_seqs = {}
		for seq in SeqIO.parse(argv[2], 'fasta'):
			if seq.id in observed:
				ref_seqs[seq.id] = str(seq.seq)

		#count coverage
		print("Subject\tReadcount\tCorrect_bases\tIncorrect_bases\tTotal_bases\tSubjlen\tCoverage\tPercent_ID")
		for o in observed:
			contig_counts = bam.count(o, start = 0, end = len(ref_seqs[o]))

			counts = bam.count_coverage(o, start = 0, end = len(ref_seqs[o]))
			#pos_ids = []
			trues = 0
			falses = 0
			total = 0
			for ref_pos in range(0, len(ref_seqs[o])):
				total += sum(counts[nt][ref_pos] for nt in range(4))

			if total == 0:
				continue

			for ref_pos in range(0, len(ref_seqs[o])):
				ref_allele = ref_seqs[o][ref_pos]
				depth = sum(counts[nt][ref_pos] for nt in range(4))
				count_a = counts[0][ref_pos]
				count_c = counts[1][ref_pos]
				count_g = counts[2][ref_pos]
				count_t = counts[3][ref_pos]
				
				if depth > 0:
					allele_counts = {'A': count_a, 'C': count_c, 'G': count_g, 'T': count_t}
					if ref_allele in allele_counts:
						false = depth - allele_counts[ref_allele]
						if false > 0:
							falses += 1
						else:
							trues += 1
					# else: N position, skip it
					
					#trues += true
					#falses += false
					#ratio = true /(true + false)
					#need the trues and positives for each ref_pos
					#print('\t'.join(str(val) for val in values) + '\t' + str(ratio))
			#pid = round(sum(pos_ids) / len(pos_ids) * 100, 2)
			#print(o + '\t' + str(contig_counts) + '\t' + str(pid))
			seqlen = len(ref_seqs[o])
			# Coverage = fraction of non-N reference positions that have any reads mapping to them. 
			# PID = fraction of non-N reference positions that have only reads matching the reference allele.
			coverage = round(((trues + falses)/seqlen) * 100, 2)
			if trues == 0 and falses == 0:
				pid = 0
			else:
				pid = round((trues/ (trues + falses)) * 100, 2)
			print(o + '\t' + str(contig_counts) + '\t'
				+ str(trues) + '\t'
				+ str(falses) + '\t'
				+ str(trues + falses) + '\t'
				+ str(seqlen) + '\t'
				+ str(coverage) + '\t'
				+ str(pid))

if __name__ == "__main__":
	main(sys.argv)
