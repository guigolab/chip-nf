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

fastqs = Channel
.from(params.input.readLines())
.map { line ->
  list = line.split()
  id = list[0]
  path = file(list[1])
  mark = list[2]
  quality = fastq(list[1]).qualityScore()
  [ id, path, mark, quality ]
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
  set prefix, file(fastq), mark, quality from fastqs

  output:
  set prefix, file("${prefix}_primary.bam"),mark into bams

  script:
  cpus = task.cpus
  memory = task.memory
  cat = fastq.name.endsWith('.gz') ? 'zcat' : 'cat'
  awk_str = 'BEGIN{OFS=FS="\\t"}$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and($2,a[i])>0){$2=xor($2,a[i])}}}{print}'
  command = ""
  command += "${cat} ${fastq} | gem-mapper -I ${index} -q offset-${quality} -T ${cpus} | pigz -p ${cpus} -c > ${prefix}.map.gz\n"
  command += "gt.filter -i ${prefix}.map.gz --max-levenshtein-error ${params.mismatches} -t ${cpus}| gt.filter --max-maps ${params.multimaps} -t ${cpus} | pigz -p ${cpus} -c > ${prefix}.filter.map.gz\n"
  command += "pigz -p ${cpus} -dc ${prefix}.filter.map.gz | gem-2-sam -T ${cpus} -I ${index} -q offset-${quality} -l --expect-single-end-reads | awk '${awk_str}' | samtools view -@ ${cpus} -Sb - | samtools sort -@ ${cpus} - ${prefix}\n"
  command += "samtools view -@ ${cpus} -bF256 ${prefix}.bam  > ${prefix}_primary.bam"
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
  set prefix, file("${prefix}_peaks*"), maxPeak into peakCallResults mode flatten
  set prefix, file("${prefix}*.bed"), maxPeak into peakCallResults mode flatten

  script:
  broad = (mark in broadMarks) ? '--broad' : ''
  command = ""
  command += "macs2 callpeak -t ${bam} -n ${prefix} --gsize hs --nomodel --extsize=${(maxPeak as int)/2} ${broad}"
  command += chipInput ? '-c ${chipInput}' : ''
}

process wiggle {

  input:
  set prefix, file(bam), maxPeak from wiggleBams

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
.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) { prefix, path, maxPeak ->
    [prefix, path, maxPeak].join("\t")
}
.subscribe {
    log.info ""
    log.info "-----------------------"
    log.info "Pipeline run completed."
    log.info "-----------------------"
}
