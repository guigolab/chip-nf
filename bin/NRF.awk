#!/usr/bin/awk -f
# compute Nonredundant fraction -  the fraction of nonredundant mapped reads in a ChIP-seq data set
{
	uniq=1
	for(i=12;i<NF;i++) {
		if ($(i)~/^XA/){
			uniq=0
			break
		}
	}
	if(uniq) {
		k=$3":"$4
		upos[k]++
		uniq++
	}
}
END {
	printf "%.4f", length(upos)/uniq
}
