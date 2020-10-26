#!/usr/bin/env python
#module load biopython/intel/python3.6/1.72
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import sys

ref = sys.argv[1]
csv = sys.argv[2]
if len(sys.argv) > 3:
	ignore_binom = sys.argv[3].lower() == 'true'
else:
	ignore_binom = False

out = csv.replace(".snplist.csv", ".vcf")

reference = list(SeqIO.parse(ref, "fasta"))[0]
ref_seq = reference.seq

chrom = "SARS-CoV2"

with open(out, 'w') as vcf:
	vcf.write("##fileformat=VCFv4.2\n")
	vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	vcf.write("##reference=file://{}\n".format(ref))
	vcf.write("##contig=<ID={},length={}>\n".format(chrom, len(ref_seq)))
	vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
	with open(csv, 'r') as infile:
		# read/skip the first line
		# saving it but don't really need it
		header_line = infile.readline()
		for line in infile:
			do_write = False
			columns = line.split(',')
			pos = int(columns[2])
			major = columns[3]
			majorfreq = float(columns[4])
			minor = columns[5]
			minorfreq = float(columns[6]) if columns[6] is not '' else ''
			totalcount = int(columns[13])
			# Binom must be True to accept a minor allele
			binom = columns[7] == "True"
			if ignore_binom: binom = True

			ref = ref_seq[pos - 1]
			
			if major == ref and minor is '': 
				continue
			elif major == "N":
				continue
			elif major == ref and minor is not '' and binom:
				alt = minor
				af = minorfreq
				dp = totalcount
				ad = "{},{}".format(round(dp*(1-af)), round(dp*af))
				do_write = True
			elif major != ref and (minor is '' or minor == ref):
				alt = major
				af = majorfreq
				dp = totalcount
				ad = "{},{}".format(round(dp*(1-af)), round(dp*af))
				do_write = True
			elif major != ref and minor is not '' and minor != ref and binom:
				alt = "{},{}".format(major, minor)
				af = "{},{}".format(majorfreq, minorfreq)
				dp = totalcount
				ad = "{},{},{}".format(round(dp*(1 - majorfreq - minorfreq)), round(majorfreq*dp), round(minorfreq*dp))
				do_write = True

			if do_write:
				info = "AF={}".format(af)
				f = "GT:AD:DP:GQ:PL	1:{}:{}:1:1".format(ad, dp)
				vcf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
					.format(chrom, pos, '.', ref, alt, '.', '.', info, f))

