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
  mark = list[3]
  quality = fastq(path).qualityScore()
  [ mergeId, id, path, mark, quality ]
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

process splitChroms {
  input:
  file genome

  output:
  file "chroms" into chromDir

  script:
  command = ""
  command += "mkdir chroms && awk '\$0~/^>/{chrom=substr(\$1,2);}{print > \"chroms/\"chrom\".fa\"}' ${genome}"
}

process mapping {
  input:
  file index from genomeIndex.val
  set mergeId, prefix, file(fastq), mark, quality from fastqs

  output:
  set mergeId, prefix, file("${prefix}_primary.bam"), mark, view into bams

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
    set mergeId, prefix, file(bam), mark, view from groupedBam

    output:
    set mergeId, prefix, file("${mergeId}.bam"), mark, view into mergedBam

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
.map { mergeId, prefix, bam, mark, view ->
  [ mergeId, bam, mark, view]
}

process model {
  input:
  set prefix, file(bam), mark, view from bams

  output:
  set prefix, file("${prefix}.params.out") into modelParams
  set prefix, file(bam), file("${prefix}.params.out"), mark, view into modelBams

  script:
  cpus = task.cpus
  command = ""
  command += "Rscript \$(which run_spp.R) -c=${bam} -rf -out=${prefix}.params.out -savp=${prefix}.pdf -p=${cpus}\n"
}

modelBams = modelBams.map { prefix, bam, paramFile, mark, view ->
  estFragLen = paramFile.text.split()[2].split(',')[0]
  [prefix, bam, mark, estFragLen, view]
}

(peakCallBams, wiggleBams, results) = modelBams.into(3)

process peakCall {
  input:
  set prefix, file(bam), mark, estFragLen, view from peakCallBams

  output:
  set prefix, file("peakOut/*"), mark, estFragLen, view into peakCallResults mode flatten

  script:
  view = "peakCall"
  broad = (mark in broadMarks) ? '--broad' : ''
  extSize = Math.round((estFragLen as int)/2)
  command = ""
  command += "macs2 callpeak -t ${bam} -n ${prefix} --outdir peakOut --gsize hs --nomodel --extsize=${extSize} ${broad}"
  command += chipInput ? '-c ${chipInput}' : ''
}

process wiggle {

  input:
  file chromSizes from chromSizes.val
  file chromDir from chromDir.val
  set prefix, file(bam), mark, estFragLen, view from wiggleBams

  output:
  set prefix, "${prefix}.bw", mark, estFragLen, view into wiggleResults

  when:
  !(params.noWiggle)

  script:
  view = 'wiggle'
  command = ""
  command += "align2rawsignal -i=${bam} -s=${chromDir} -u=${genomeMapDir}"
  command += " -o=${prefix}.bedgraph -of=bg -v=stdout -l=${(estFragLen as int)-50}"
  command += " -mm=${Math.round(task.memory.toGiga() * 0.5)}\n"
  command += "bedGraphToBigWig ${prefix}.bedgraph ${chromSizes} ${prefix}.bw"
}

results.mix(peakCallResults, wiggleResults)
.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) { prefix, path, mark, maxPeak, view ->
    [ prefix, path, mark, maxPeak, view ].join("\t")
}
.subscribe {
    log.info ""
    log.info "-----------------------"
    log.info "Pipeline run completed."
    log.info "-----------------------"
}
