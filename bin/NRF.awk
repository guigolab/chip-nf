#!/usr/bin/awk -f
# compute Nonredundant fraction -  the fraction of nonredundant mapped reads in a ChIP-seq data set
{
	for(i=12;i<NF;i++) {
		if ($(i)=="NH:i:1"){
			k=$3":"$4
			upos[k]++
			uniq++
			break
		}
	}
}
END {
	printf "%.4f", length(upos)/uniq
}