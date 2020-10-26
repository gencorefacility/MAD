#!/usr/bin/env python3

import sys

# give it a vcf file
in_vcf = sys.argv[1]

# give it a DP
dp = int(float(sys.argv[2]))

with open(in_vcf, 'r') as vcf:
	for line in vcf:
		if line.startswith('#'):
			print(line.strip())
		else:
			columns = line.split('\t')
			pos = columns[1]
			ref = columns[3]
			alt = columns[4]
			chrom = columns[0]
			af = float(columns[7].split("=")[1]	)
			format_field = columns[8]

			sample_field = columns[9]
			gt = sample_field.split(':')[0]

			new_ad_1 = round(dp*af)
			new_ad_0 = dp - new_ad_1
			new_sample_field = "{}:{},{}:{}:1:1".format(gt, new_ad_0, new_ad_1, dp)			

			print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
				.format(chrom, pos, '.', ref, alt, '.', '.', columns[7], format_field, new_sample_field))

