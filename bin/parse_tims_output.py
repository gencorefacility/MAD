#!/usr/bin/env python3
import sys
from pyfaidx import Fasta
import glob

fa = sys.argv[1]
sample_id = sys.argv[2]
if len(sys.argv) > 3:
	ignore_binom = sys.argv[3].lower() == 'true'
else:
	ignore_binom = False

contigs = Fasta(fa)

out = "{}.vcf".format(sample_id)

with open(out, 'w') as vcf:
	vcf.write("##fileformat=VCFv4.2\n")
	vcf.write("##source=timo\n")
	vcf.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">\n')
	vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	vcf.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Number of observation for each allele">\n')
	vcf.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
	vcf.write('##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype">\n')
	vcf.write("##reference=file://{}\n".format(fa))

	vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
	for chrom in contigs.keys():
		csv = glob.glob("FILES/fullvarlist/filtered.STRAIN.{}.*.snplist.csv".format(chrom))[0]
		ref_seq = contigs[chrom]
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

