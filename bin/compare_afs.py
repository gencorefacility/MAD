#!/usr/bin/env python3
# requires: pysam/intel/python3.6/0.14.1
# Expects you to have run the command
# `bcftools isec $golden_vcf $workflow_vcf -p isec`
# Expects 3 files 000[123].vcf in a dir called isec

import pysam
import sys
import json
import csv
import argparse

def main():
	## Retrieve command line arguments
	in_arg = get_input_args()
	golden_vcf = in_arg.golden_vcf
	workflow_vcf = in_arg.workflow_vcf

	# Output from bcftools isec
	false_negatives_vfc = 'isec/0000.vcf'
	false_positives_vcf = 'isec/0001.vcf'
	true_positives_vcf = 'isec/0002.vcf'

	# These are lists which hold the positions of the fn, fp, and tp vars
	false_negatives = get_positions_from_vcf(false_negatives_vfc)
	false_positives = get_positions_from_vcf(false_positives_vcf)
	true_positives = get_positions_from_vcf(true_positives_vcf)

	results = []

	vcf_golden = pysam.VariantFile(golden_vcf)
	vcf_workflow = pysam.VariantFile(workflow_vcf)
	
	sample_id = workflow_vcf.split(".vcf")[0]

	for pos in false_negatives:
		data = get_data_golden(vcf_golden, pos)
		#print(data)
		results.append({
			'sample_id': sample_id, 
			'pos': pos, 
			'af_golden': data['af'], 
			'af_workflow': 0, 
			'ref': data['ref'], 
			'alt': data['alt'], 
			'dp': data['dp']
		})

	for pos in false_positives:
		data = get_data_workflow(vcf_workflow, pos)
		results.append({
			'sample_id': sample_id, 
			'pos': pos, 
			'af_golden': 0, 
			'af_workflow': data['af'], 
			'ref': data['ref'], 
			'alt': data['alt'], 
			'dp': data['dp']
		})

	for pos in true_positives:
		gold_data = get_data_golden(vcf_golden, pos)
		workflow_data = get_data_workflow(vcf_workflow, pos)
		results.append({
			'sample_id': sample_id, 
			'pos': pos, 
			'af_golden': gold_data['af'], 
			'af_workflow': workflow_data['af'], 
			'ref': workflow_data['ref'], 
			'alt': workflow_data['alt'], 
			'dp': workflow_data['dp']
		})

	with open("{}.csv".format(in_arg.out), 'w', newline='') as f:
		w = csv.DictWriter(f, fieldnames = results[0].keys())
		w.writeheader()
		w.writerows(results)

def get_input_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('--golden_vcf', type = str, required = True)
	parser.add_argument('--workflow_vcf', type = str, required = True)
	parser.add_argument('--out', type = str, required = True)
	in_arg = parser.parse_args()
	return in_arg

def get_positions_from_vcf(vcf):
	l = []
	vcf = pysam.VariantFile(vcf, "r")
	for var in vcf:	
		l.append(var.pos)
	return l

def get_data_workflow(vcf, pos):
	# Parse data from variant caller VCFs
	# assumes only 1 variant at this position, takes first one
	var = list(vcf.fetch(vcf.get_reference_name(0), pos - 1, pos))[0]
	# Set dp = None as default
	dp = None
	for sample in var.samples:
		ad = var.samples[sample]['AD']
		dp = int(var.samples[sample]['DP'])

	# Get AF directly from INFO AF field for some tools
	use_af_tool_list = ['varscan', 'ivar', 'tim', 'cliquesnv.vcf', 'lofreq.vcf']
	if any(x in vcf.filename.decode() for x in use_af_tool_list):
		info_af = var.info["AF"]
		if type(info_af) is tuple:
			af_workflow = round(float(var.info["AF"][0]),2)
		elif type(info_af) in [float, str]:
			af_workflow = round(float(var.info["AF"]),2)
	else:
		for sample in var.samples:
	                ad = var.samples[sample]['AD']
        	        dp = int(var.samples[sample]['DP'])
		af_workflow = int(ad[1]) / dp
	
	return {"af": af_workflow, "dp": dp, "ref": var.ref, "alt": var.alts[0]}

def get_data_golden(vcf, pos):
        var = list(vcf.fetch(vcf.get_reference_name(0), pos - 1, pos))[0]
        for sample in var.samples:
                dp = int(var.samples[sample]['DP'])
        af_golden = round(float(var.info["AF"][0]),2)
        return {"af": af_golden, "dp": dp, "ref": var.ref, "alt": var.alts[0]}

if __name__ == "__main__":
	main()

