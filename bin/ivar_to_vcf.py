#!/usr/bin/env python
import sys

vs_out = sys.argv[1]

out = vs_out.replace(".tsv", ".vcf")

if len(sys.argv) > 2:
        include_non_pf = sys.argv[3].lower() == 'true'
else:
        include_non_pf = False

with open(out, 'w') as vcf:
	vcf.write("##fileformat=VCFv4.2\n")
	vcf.write("##source=iVar\n")
	vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
	vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n")
	with open(vs_out, 'r') as infile:
		# read/skip the first line
		# saving it but don't really need it
		header_line = infile.readline()
		for line in infile:
			columns = line.split('\t')
			passed = columns[13] == "TRUE"
			if not passed and not include_non_pf:
				continue
			chrom = columns[0]
			pos = int(columns[1])
			ref = columns[2]
			af = float(columns[10])
			alt = columns[3]
			if len(alt) > 1: 
				continue
			
			dp = int(columns[11])
			ad = "{},{}".format(round(dp*(1-af)), round(dp*af))

			info = "AF={}".format(af)
			f = "GT:AD:DP:GQ:PL	1:{}:{}:1:1".format(ad, dp)
			
			vcf.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
					.format(chrom, pos, '.', ref, alt, '.', '.', info, f))
