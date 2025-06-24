# usage: python3 post_who_gc_content.py <whokaryote_predictions.tsv> <assembly.fasta> <out_contigs.tsv>
# Can also be used as Snakemake script

import sys
import gzip
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def parse_tsv(input_tsv):
	contig_dict = {}
	with open(input_tsv, 'r') as infile:
		infile.readline()
		for line in infile:
			contig, who = line.strip().split('\t')
			contig_dict[contig] = {
				'who': who,
				'gc': 0
			}
	return contig_dict


def parse_fasta(input_fasta, contig_dict):
	with open(input_fasta, 'r') as infile:
		for record in SeqIO.parse(infile, 'fasta'):
			try:
				contig_dict[record.id]['gc'] = gc_fraction(record.seq) * 100
			except KeyError:
				print(record.id, len(record.seq))
	return contig_dict

def parse_fastq(input_fastq):
	reads_dict = {}
	with gzip.open(input_fastq, 'rt') as infile:
		for record in SeqIO.parse(infile, 'fastq'):
			reads_dict[record.id] = gc_fraction(record.seq) * 100

	return reads_dict


if __name__ == '__main__':
	# Check if running as Snakemake script
	if 'snakemake' in globals():
		# Get input/output from Snakemake
		whokaryote_tsv = snakemake.input[0]  # whokaryote predictions
		assembly_fasta = snakemake.input[1]  # assembly fasta
		output_tsv = snakemake.output[0]     # output contigs with GC
		
		# Parse data
		contig_dict = parse_tsv(whokaryote_tsv)
		contig_dict = parse_fasta(assembly_fasta, contig_dict)
		
		# Write output
		with open(output_tsv, 'w') as outfile:
			outfile.write("contig\twho\tgc_content\n")  # Header
			for contig in contig_dict:
				outfile.write(f'{contig}\t{contig_dict[contig]["who"]}\t{contig_dict[contig]["gc"]}\n')
	
	else:
		# Command line execution (fallback)
		if len(sys.argv) < 4:
			print("Usage: python3 post_who_gc_content.py <whokaryote_predictions.tsv> <assembly.fasta> <out_contigs.tsv> [raw_reads.fastq.gz] [out_reads.tsv]")
			sys.exit(1)
		
		contig_dict = parse_tsv(sys.argv[1])
		contig_dict = parse_fasta(sys.argv[2], contig_dict)
		
		# Write contigs output
		with open(sys.argv[3], 'w') as outfile:
			outfile.write("contig\twho\tgc_content\n")  # Header
			for contig in contig_dict:
				outfile.write(f'{contig}\t{contig_dict[contig]["who"]}\t{contig_dict[contig]["gc"]}\n')
		
		# Optional reads processing (if provided)
		if len(sys.argv) >= 6:
			reads_dict = parse_fastq(sys.argv[4])
			with open(sys.argv[5], 'w') as outfile:
				outfile.write("read\tgc_content\n")  # Header
				for read in reads_dict:
					outfile.write(f'{read}\t{reads_dict[read]}\n')