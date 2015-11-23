params.mismatches = 2
params.multimaps = 10
params.dbFile = 'chipseq-pipeline.db'

chipInput = null

broadMarks = [
  "H3K27me3",
  "H3K36me3",
  "H3K9me3",
  "H3K4me1"
]

pdb = file(params.dbFile)
pdb.write('')

genome = file(params.genome)
input = file(params.input)
genomeMapDir = file(params.genomeMapDir)

fastqs = Channel
.from(input.readLines())
.map { line ->
  list = line.split()
  mergeId = list[0]
  id = list[1]
  path = file(list[2])
  controlId = list[3]
  mark = list[4]
  quality = fastq(path).qualityScore()
  [ mergeId, id, path, controlId, mark, quality ]
}

process index {
  input:
  file genome

  output:
  file "genome_index.gem" into genomeIndex

  script:
  command = ""
  command += "sed 's/ .*//' ${genome} > genome_processed.fa\n"
  command += "gem-indexer -i genome_processed.fa -o genome_index -T ${task.cpus} -m ${task.memory.toBytes()} && rm genome_processed.fa"
}

process fastaIndex {
  input:
  file genome

  output:
  file "${genome}.fai" into chromSizes

  script:
  command = ""
  command += "samtools faidx ${genome}"
}

process mapping {
  input:
  file index from genomeIndex.val
  set mergeId, prefix, file(fastq), controlId, mark, quality from fastqs

  output:
  set mergeId, prefix, file("${prefix}_primary.bam"), controlId, mark, view into bams

  script:
  cpus = task.cpus
  memory = task.memory
  readGroup = "ID=${prefix},SM=${mergeId}"
  view = "Alignments"
  cat = fastq.name.endsWith('.gz') ? 'zcat' : 'cat'
  awk_str = 'BEGIN{OFS=FS="\\t"}$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and($2,a[i])>0){$2=xor($2,a[i])}}}{print}'
  command = ""
  command += "${cat} ${fastq} | gem-mapper -I ${index} -q offset-${quality} -T ${cpus} | pigz -p ${cpus} -c > ${prefix}.map.gz\n"
  command += "gt.filter -i ${prefix}.map.gz --max-levenshtein-error ${params.mismatches} -t ${cpus}| gt.filter --max-maps ${params.multimaps} -t ${cpus} | pigz -p ${cpus} -c > ${prefix}.filter.map.gz\n"
  command += "pigz -p ${cpus} -dc ${prefix}.filter.map.gz | gem-2-sam -T ${cpus} -I ${index} -q offset-${quality} -l --expect-single-end-reads --read-group ${readGroup} | awk '${awk_str}' | samtools view -@ ${cpus} -Sb - | samtools sort -@ ${cpus} - ${prefix}\n"
  command += "samtools view -@ ${cpus} -bF256 ${prefix}.bam  > ${prefix}_primary.bam"
}

// Merge or rename bam
singleBam = Channel.create()
groupedBam = Channel.create()

bams.groupTuple(by: [0,3,4])
.choice(singleBam, groupedBam) {
  it[2].size() > 1 ? 1 : 0
}

process mergeBam {

    input:
    set mergeId, prefix, file(bam), controlId, mark, view from groupedBam

    output:
    set mergeId, prefix, file("${mergeId}.bam"), controlId, mark, view into mergedBam

    script:
    cpus = task.cpus
    prefix = prefix.sort().join(':')
    command = ""
    command += """(
      samtools view -H ${bam} | grep -v "@RG";
	    for f in ${bam}; do
		    samtools view -H \$f | grep "@RG";
	    done) > header.txt \n"""
    command += "samtools merge -@ ${cpus} -h header.txt ${mergeId}.bam ${bam}"

}


bams = singleBam
.mix(mergedBam)
.map { mergeId, prefix, bam, controlId, mark, view ->
  [ mergeId, bam, controlId, mark, view]
}

// cross bams and inputs
treat = Channel.create()
control = Channel.create()
bams.choice(treat, control) {
    it[3] == 'input' ? 1 : 0
}

bams = control.cross(treat) { it[2] }.map { c, t ->
    [t[0], t[1], c[1], t[3], t[4]]
}

process model {
  input:
  set prefix, file(bam), file(control), mark, view from bams

  output:
  set prefix, file("${prefix}.params.out") into modelParams
  set prefix, file(bam), file(control), file("${prefix}.params.out"), mark, view into modelBams

  script:
  cpus = task.cpus
  command = ""
  command += "Rscript \$(which run_spp.R) -c=${bam} -rf -out=${prefix}.params.out -savp=${prefix}.pdf -p=${cpus}\n"
}

modelBams = modelBams.map { prefix, bam, control, paramFile, mark, view ->
  fragLen = paramFile.text.split()[2].split(',')[0]
  [prefix, bam, control, mark, fragLen, view]
}

(peakCallBams, results) = modelBams.into(2)

process peakCall {
  input:
  file chromSizes
  set prefix, file(bam), file(control), mark, fragLen, view from peakCallBams

  output:
  set prefix, file("peakOut/${prefix}_peaks.narrowPeak"), mark, estFragLen, val("narrowPeak") into peakCallResults
  set prefix, file("peakOut/${prefix}_peaks.broadPeak"), mark, estFragLen, val("broadPeak") into peakCallResults
  set prefix, file("peakOut/${prefix}_peaks.gappedPeak"), mark, estFragLen, val("gappedPeak") into peakCallResults
  set prefix, file("peakOut/${prefix}.fc_signal.bw"), mark, estFragLen, val("fcSignal") into peakCallResults
  set prefix, file("peakOut/${prefix}.pval_signal.bw"), mark, estFragLen, val("pvalueSignal") into peakCallResults

  script:
  //extSize = Math.round((fragLen as int)/2)
  command = ""
  // narrow peaks and preliminary signal tracks
  command += "macs2 callpeak -t ${bam} -c ${control} -n ${prefix} --outdir peakOut"
  command += " -f BAM -g hs -p 1e-2 --nomodel --extsize=${fragLen}"
  command += " --keep-dup all -B --SPMR\n"
  // Broad and Gapped Peaks
  command += "macs2 callpeak -t ${bam} -c ${control} -n ${prefix} --outdir peakOut"
  command += " -f BAM -g hs -p 1e-2 --broad --nomodel --extsize=${fragLen}"
  command += " --keep-dup all\n"
  // rescale peaks on 10-1000 scale
  if ( params.rescale ) {
    ['narrow', 'broad', 'gapped'].collect { type ->
      rescale_awk_str = 'BEGIN{FS=OFS="\\t";min=1e20;max=-1}'
      rescale_awk_str += 'NR==FNR&&NF!=0{min>$5?min=$5:min=min;max<$5?max=$5:max=max;next}'
      rescale_awk_str += 'NF!=0{n=$5;x=10;y=1000;$5=int(((n-min)*(y-x)/(max-min))+x);print}'
      command += "awk '${rescale_awk_str}' peakOut/${prefix}_peaks.${type}Peak peakOut/${prefix}_peaks.${type}Peak"
      command += " > peakOut/${prefix}_peaks.${type}Peak_rescaled && mv peakOut/${prefix}_peaks.${type}Peak{_rescaled,}\n"
    }
  }
  // Fold enrichment signal tracks
  command += "macs2 bdgcmp -t peakOut/${prefix}_treat_pileup.bdg"
	command += " -c peakOut/${prefix}_control_lambda.bdg --outdir peakOut"
  command += " -o ${prefix}_FE.bdg -m FE\n"
  command += "slopBed -i peakOut/${prefix}_FE.bdg -g ${chromSizes} -b 0"
  command += " | bedClip stdin ${chromSizes} peakOut/${prefix}.fc.signal.bedgraph\n"
  command += "bedGraphToBigWig peakOut/${prefix}.fc.signal.bedgraph ${chromSizes} peakOut/${prefix}.fc_signal.bw\n"
  // -log10(p-value) signal tracks
  command += "bamReads=\$(samtools view -c ${bam}) && controlReads=\$(samtools view -c ${control})"
  command += " && sval=\$(bc <<< \$((bamReads<controlReads?bamReads:controlReads))/1000000)\n"
  command += "macs2 bdgcmp -t peakOut/${prefix}_treat_pileup.bdg"
	command += " -c peakOut/${prefix}_control_lambda.bdg --outdir peakOut"
  command += " -o ${prefix}_ppois.bdg -m ppois -S \$sval\n"
  command += "slopBed -i peakOut/${prefix}_ppois.bdg -g ${chromSizes} -b 0"
  command += " | bedClip stdin ${chromSizes} peakOut/${prefix}.pval.signal.bedgraph\n"
  command += "bedGraphToBigWig peakOut/${prefix}.pval.signal.bedgraph ${chromSizes} peakOut/${prefix}.pval_signal.bw\n"
  command += "rm -rf peakOut/${prefix}*.bdg peakOut/${prefix}*.bedgraph peakOut/${prefix}*.xls peakOut/${prefix}*.bed"
}

results.map { prefix, bam, control, mark, fragLen, view ->
  [ prefix, bam, mark, fragLen, view ]
}
.mix(peakCallResults)
.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) { prefix, path, mark, fragLen, view ->
    [ prefix, path, mark, fragLen, view ].join("\t")
}
.subscribe {
    log.info ""
    log.info "-----------------------"
    log.info "Pipeline run completed."
    log.info "-----------------------"
}
