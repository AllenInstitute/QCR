#!/bin/bash
data_path=$1"/*labelled.RDS"
file_names=`ls $data_path`

echo -e "samples:"	
echo -e ' \t '; 
for i in $file_names
do
	j="${i%.*}"
	k=$(basename $j)
	sample=$(echo $k | cut -d"_" -f1)
	echo -e '\t' ${sample} ":" ${sample}
done
echo "workdir:"
echo "reference:"
echo "resolution:"
echo "ref_subclass:"
echo "ref_class:"
echo "sample_id:"
echo "donor:"
echo "exclude:"
echo "glia_cutoff:"
echo "neuron_cutoff:"
echo "ref_subclass_color:"
