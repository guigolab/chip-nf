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

// fastqPattern = ~"(.+)_[0-9]+_[ACGTN]+.fastq(.gz)?"

chromSizes="/users/rg/projects/references/Genome/H.sapiens/hg19/hg19.chromsizes"
chromDir="/users/rg/projects/ERC/chipseq-human/genome/chrs"
genomeMapDir="/users/rg/projects/references/Genome/H.sapiens/hg19/globalmap_k20tok54"

index = file(params.index)
input = file(params.input)

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
// fastqs = Channel
// .fromPath(params.input).map {
//     matcher = it.name =~ fastqPattern
//     if (!matcher) {
//       exit 1, "Cannot match fastq name"
//     }
//     [matcher.group(1), it, fastq(it).qualityScore()]
// }

/*process index {
  input:
  file genome

  """
  gem-indexer -i {genome} -o ${genome_pref} -T ${cpus} -m ${memory}
  """
}*/

process mapping {
  input:
  file index
  set mergeId, prefix, file(fastq), mark, quality from fastqs

  output:
  set mergeId, prefix, file("${prefix}_primary.bam"), mark into bams

  script:
  cpus = task.cpus
  memory = task.memory
  readGroup = "ID=${prefix},SM=${mergeId}"
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

bams.groupTuple(by: [0,3])
.choice(singleBam, groupedBam) {
  it[2].size() > 1 ? 1 : 0
}

process mergeBam {

    input:
    set mergeId, prefix, file(bam), mark from groupedBam

    output:
    set mergeId, prefix, file("${mergeId}.bam"), mark into mergedBam

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
.map { mergeId, prefix, bam, mark ->
  [ mergeId, bam, mark ]
}

process model {
  input:
  set prefix, file(bam), mark from bams

  output:
  set prefix, file("${prefix}.params.out") into modelParams
  set prefix, file(bam), file("${prefix}.params.out"), mark into modelBams

  script:
  cpus = task.cpus
  command = ""
  command += "Rscript \$(which run_spp.R) -c=${bam} -rf -out=${prefix}.params.out -savp=${prefix}.pdf -p=${cpus}\n"
}

modelBams = modelBams.map { prefix, bam, paramFile, mark ->
  maxPeak = paramFile.text.split()[2].split(',')[0]
  [prefix, bam, mark, maxPeak]
}

(peakCallBams, wiggleBams, results) = modelBams.into(3)

process peakCall {
  input:
  set prefix, file(bam), mark, maxPeak from peakCallBams

  output:
  set prefix, file("peakOut/*"), mark, maxPeak into peakCallResults mode flatten

  script:
  broad = (mark in broadMarks) ? '--broad' : ''
  command = ""
  command += "macs2 callpeak -t ${bam} -n ${prefix} --gsize hs --nomodel --extsize=${(maxPeak as int)/2} --outdir peakOut ${broad}"
  command += chipInput ? '-c ${chipInput}' : ''
}

process wiggle {

  input:
  set prefix, file(bam), mark, maxPeak from wiggleBams

  output:
  set prefix, "${prefix}.bw", maxPeak into wiggleResults

  when:
  params.wiggle

  script:
  command = ""
  command += "align2rawsignal -i=${bam} -s=${chromDir} -u=${genomeMapDir} -o=${prefix}.bedgraph -of=bg -v=stdout -l=${(maxPeak as int)-50}\n"
  command += "bedGraphToBigWig ${prefix}.bedgraph ${chromSizes} ${prefix}.bw"
}

results.mix(peakCallResults, wiggleResults)
.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) { prefix, path, mark, maxPeak ->
    [ prefix, path, mark, maxPeak ].join("\t")
}
.subscribe {
    log.info ""
    log.info "-----------------------"
    log.info "Pipeline run completed."
    log.info "-----------------------"
}
