#!/bin/bash
# TODO:: fix this
# This file is used to convert mulitiple RGI json files into tab-delimited outputs
# usage: bash rgi_jsontab_batch.bash <DIRECTORY>

echo $1

if [ -d "$1" ]; then
  for f in $1; do
  	echo ${f##*/}
	g=${f##*/};
	echo rgi_jsontab -i ${f##*/} -o ${g%.*};
done
else
	echo "path: '$1' does not exist"
fi