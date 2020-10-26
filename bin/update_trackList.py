#!/usr/bin/env python3

import json
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('--fcid', type = str, required = True)
parser.add_argument('--id', type = str, required = True)
parser.add_argument('--bw', action='store_true', default=False)
parser.add_argument('--regex', type = str, required = False)
in_arg = parser.parse_args()

sample_id = in_arg.id
fcid = in_arg.fcid

# If we specify a regex, we use this to clean up 
# the sample_id for the category
if in_arg.regex is not None:
	category = re.sub(in_arg.regex, "", sample_id)
else:
	category = sample_id

snps = {
	"category": category,
	"key": sample_id + "-snps",
	"label": sample_id + "-snps",
	"storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
	"type": "CanvasFeatures",
	"urlTemplate": "tracks/samples/{}/{}/{}_filtered_snps.ann.vcf.gz".format(fcid, category, sample_id),
	"style": {
		"color": "{variantColor}"
	}
}

indels = {
        "category": category,
        "key": sample_id + "-indels",
        "label": sample_id + "-indels",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        "type": "CanvasFeatures",
        "urlTemplate": "tracks/samples/{}/{}/{}_filtered_indels.vcf.gz".format(fcid, category, sample_id),
	"style": {
		"color": "{variantColor}"
	}
}

bam = {
        "category": category,
        "key": sample_id + "-bam",
        "label": sample_id + "-bam",
	"chunkSizeLimit": 40000000,
        "storeClass": "JBrowse/Store/SeqFeature/BAM",
        "type": "JBrowse/View/Track/Alignments2",
        "urlTemplate": "tracks/samples/{}/{}/{}_downsampled.bam".format(fcid, category, sample_id)
}

bw = {
        "category": category,
        "key": sample_id + "-coverage",
        "label": sample_id + "-coverage",
        "storeClass": "JBrowse/Store/SeqFeature/BigWig",
        "type": "JBrowse/View/Track/Wiggle/XYPlot",
        "urlTemplate": "tracks/samples/{}/{}/{}_coverage.bam.bw".format(fcid, category, sample_id),
	"autoscale" : "local",
	"scale" : "log",
	"bg_color" : "white",
	"variance_band" : True,
	"style" : {
		"height" : 50
	}
}

snp_cov = {
	"category": category,
	"key": sample_id + "-snpcoverage-freebayes",
	"label": sample_id + "-snpcoverage-freebayes",
	"storeClass": "MPileupPlugin/Store/SeqFeature/BigWigVcfAdDp",
	"type": "MPileupPlugin/View/Track/SNPCoverage",
	"bigwig": {
		"storeClass": "JBrowse/Store/SeqFeature/BigWig",
		"urlTemplate": "tracks/samples/{}/{}/{}_coverage.bam.bw".format(fcid, category, sample_id)
	},
	"vcf": {
		"storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
		"urlTemplate": "tracks/samples/{}/{}/{}_freebayes.vcf.gz".format(fcid, category, sample_id)
	}
}

ivar_snps = {
        "category": category,
        "key": sample_id + "-ivar",
        "label": sample_id + "-ivar",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        "type": "CanvasFeatures",
        "urlTemplate": "tracks/samples/{}/{}/{}_ivar.vcf.gz".format(fcid, category, sample_id),
	"style": {
		"color": "{variantColor}"
	}
}

timo_snps = {
        "category": category,
        "key": sample_id + "-timo-snps",
        "label": sample_id + "-timo-snps",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        "type": "CanvasFeatures",
        "urlTemplate": "tracks/samples/{}/{}/{}_tims.vcf.gz".format(fcid, category, sample_id),
	"style": {
		"color": "{variantColor}"
	}
}

hc_snps = {
        "category": category,
        "key": sample_id + "-hc-snps",
        "label": sample_id + "-hc-snps",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        "type": "CanvasFeatures",
        "urlTemplate": "tracks/samples/{}/{}/{}_filtered_snps.vcf.gz".format(fcid, category, sample_id),
        "style": {
                "color": "{variantColor}"
        }
}

mutect2_snps = {
        "category": category,
        "key": sample_id + "-mutect2-snps",
        "label": sample_id + "-mutect2-snps",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        "type": "CanvasFeatures",
        "urlTemplate": "tracks/samples/{}/{}/{}_mutect2_filtered.vcf.gz".format(fcid, category, sample_id),
        "style": {
                "color": "{variantColor}"
        }
}

freebayes_default_snps = {
        "category": category,
        "key": sample_id + "-freebayes-default-snps",
        "label": sample_id + "-freebayes-default-snps",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        "type": "CanvasFeatures",
        "urlTemplate": "tracks/samples/{}/{}/{}_freebayes_default.vcf.gz".format(fcid, category, sample_id),
        "style": {
                "color": "{variantColor}"
        }
}

varscan_snps = {
        "category": category,
        "key": sample_id + "-varscan-snps",
        "label": sample_id + "-varscan-snps",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        "type": "CanvasFeatures",
        "urlTemplate": "tracks/samples/{}/{}/{}_varscan.vcf.gz".format(fcid, category, sample_id),
        "style": {
                "color": "{variantColor}"
        }
}

golden_snps = {
        "category": category,
        "key": sample_id + "-golden-snps",
        "label": sample_id + "-golden-snps",
        "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
        "type": "CanvasFeatures",
        "urlTemplate": "tracks/samples/{}/{}/{}_golden.vcf.gz".format(fcid, category, sample_id),
        "style": {
                "color": "{variantColor}"
        }
}

new_trackList = "{}_trackList.json".format(fcid)
with open(new_trackList) as json_file:
	data = json.load(json_file)

tracks = data['tracks']
	
if not in_arg.bw:
	tracks.append(bam)
	tracks.append(snps)
	tracks.append(ivar_snps)
	tracks.append(timo_snps)
	tracks.append(snp_cov)
	#tracks.append(hc_snps)
	tracks.append(mutect2_snps)
	tracks.append(freebayes_default_snps)
	tracks.append(varscan_snps)
	tracks.append(golden_snps)
tracks.append(bw)

with open(new_trackList, 'w') as f: 
	json.dump(data, f, indent=4) 
