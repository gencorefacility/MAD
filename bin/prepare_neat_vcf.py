#!/usr/bin/env python3

import sys
import random
import numpy as np
import pysam
import os
from pyfaidx import Fasta

# give it a vcf file
in_vcf = sys.argv[1]

# give it a ref fa

fa = sys.argv[2]

# give it an AF
in_af = sys.argv[3]

vcf = pysam.VariantFile(in_vcf, "r")
assembly = os.path.basename(fa)
contigs = Fasta(fa)

def get_af(af):
	if af == 'random':
		vals = np.arange(0.01, 1.01, 0.01).tolist()
		p_under_50 = 3
		p_over_50 = 1
		l1 = [p_under_50] * 50
		l2 = [p_over_50] * 50	
		weights = l1 + l2
		random_af = random.choices(vals, weights=weights, k=1)[0]
		return random_af
	else:
		af = float(af)
		if af > 1:
			print("AF cannot be greater than 1! AF = ", af)
			sys.exit()
		return float(af)

contigs_added = False	
with open(in_vcf, 'r') as vcf:
	for line in vcf:
		if line.startswith('#'):
			if line.startswith("#CHROM\tPOS"):
				print(line.strip() + "\tFORMAT\tsample")
			else:
				print(line.strip())

			## Add Contigs Here (needed by ReSeq)
			if not contigs_added:
				for key in contigs.keys():
					print("##contig=<ID={},length={},assembly={}>".format(key,len(contigs[key]),assembly))
				contigs_added = True
		else:
			columns = line.split('\t')
			pos = columns[1]
			ref = columns[3]
			alt = columns[4]
			chrom = columns[0]

			af = get_af(in_af)
			
			new_gt_field = 'GT:AD:DP:GQ:PL\t1'
			# - 1 in the range function because we've
			# already added '1' in the line above
			for i in range(int(af*float(100)) - 1):
				new_gt_field+= '/1'
			# compare_afs.py will use ad[1]/dp to calculate the AF
			# set it appropriately here. Using dp = 100 for simplicity
			dp = 100
			new_ad_1 = round(dp*af)
			new_ad_0 = dp - new_ad_1
			new_gt_field+= ':{},{}:{}:1:1'.format(new_ad_0, new_ad_1, dp)
			print("{}\t{}\t{}\t{}\t{}\t{}\t{}\tAF={}\t{}"
				.format(chrom, pos, '.', ref, alt, '.', '.', af, new_gt_field))

