#!/bin/bash
pair_ids=$1
fcid=$2

echo "pair_ids: $pair_ids"
echo "fcid: $fcid"

# Create a new trackList for this fcid
# update_trackList.py called below will look
# for this file and append tracks to it
echo '{"formatVersion": 1, "tracks": []}' > ${fcid}_trackList.json

# For pairs (ex: A + B, unmerged), just update the track info
# Data will be sent altogether into the pair_id dir (next)
for pair_id in ${pair_ids}
do
	# pair_id might have a comma, or leading
	# or trailing square bracket, get rid of it
	pair_id=$(echo $pair_id | sed "s/[][,]//g")

	# create the dir where we will send all the data for this
	# pair_id then rsync all the data for this pair_id there
	dir=/usr/share/nginx/html/dev/data/tracks/samples/${fcid}/${pair_id}
	ssh jbrowse@covid-19.bio.nyu.edu "mkdir -p ${dir}"
	rsync --copy-links ${pair_id}[-_]* \
		jbrowse@covid-19.bio.nyu.edu:${dir}/.

	# update the trackList for all pair_id (merged) files (default)
	update_trackList.py --fcid="$fcid" --id="$pair_id"
done

