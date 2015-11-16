params.mismatches = 2
params.multimaps = 10

broadPeaks = [
  "H3K27me3",
  "H3K36me3",
  "H3K9me3",
  "H3K4me1"
]

fastqPattern = ~"(.+)_[0-9]+_[ACGTN]+.fastq(.gz)?"

chromSizes="/users/rg/projects/references/Genome/H.sapiens/hg19/hg19.chromsizes"
chromDir="/users/rg/projects/ERC/chipseq-human/genome/chrs"
genomeMapDir="/users/rg/projects/references/Genome/H.sapiens/hg19/globalmap_k20tok54"

index = file(params.index)
fastqs = Channel
.fromPath(params.input).map {
    matcher = it.name =~ fastqPattern
    if (!matcher) {
      exit 1, "Cannot match fastq name"
    }
    [matcher.group(1), it, fastq(it).qualityScore()]
}

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
  set prefix, file(fastq), quality from fastqs

  output:
  set prefix, file("${prefix}_primary.bam") into bams

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

(modelBams, wiggleBams, peakCallBams, results) = bams.into(4)

process model {
  input:
  set prefix, bam from modelBams

  output:
  // stdout in peaks
  set prefix, file("${prefix}.params.out") into modelParams

  script:
  cpus = task.cpus
  prefix = bam.baseName
  command = ""
  command += "Rscript \$(which run_spp.R) -c=${bam} -rf -out=${prefix}.params.out -savp=${prefix}.pdf -p=${cpus}\n"
}

modelParams = modelParams.map { prefix, out ->
  maxPeak = out.text.split()[2].split(',')[0]
  [prefix, out, maxPeak]
}

process peakCall {
  input:
  set prefix, file(bam)bam from peakCallBams

  output:
  set file("${prefix}_peaks.xls"), file("${prefix}_summits.bed"), file("${prefix}_peaks.*Peak") into results

  script:
  broad = (peak in broad_peaks) ? '--broad' : ''
  command = ""
  command += "macs2 callpeak -t ${bam} -n ${prefix} --gsize hs -c ${chip_input} --nomodel --shiftsize=half_fragment_size ${broad}"
}

process wiggle {

  input:
  set prefix, bam from wiggleBams
  set prefix, out, maxPeak from modelParams

  output:
  set prefix, "${prefix}.bw", maxPeak into results

  script:
  command = ""
  command += "align2rawsignal -i=${bam} -s=${chromDir} -u=${genomeMapDir} -o=${prefix}.bedgraph -of=bg -v=stdout -l=${maxPeak-50}\n"
  command += "bedGraphToBigWig ${prefix}.bedgraph ${chromSizes} ${prefix}.bw"
}

results.println { it }
