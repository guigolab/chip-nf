params.mismatches = 2
params.multimaps = 10

index = params.index
fastqs = Channel.from(params.input)

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
  file fastq from fastqs

  output:
  file ${out_prefix}_primary.bam into bams

  script:
  cat = { fastq.name.endsWith('.gz') ? 'zcat' : 'cat' }
  fqpattern = ~"(.+)_[0-9]+_[ACGTN]+.fastq(.gz)?"
  matcher = fastq.name =~ fqpattern
  if (!matcher) {
    exit 1, "Cannot match fastq name"
  }
  out_prefix = matcher.group(1)
  awk_str = 'BEGIN{OFS=FS=\"\\t\"}$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and($2,a[i])>0){$2=xor($2,a[i])}}}{print}\'")'
  command = ""
  command += "${cat} ${fastq} | gem-mapper -I ${index} -q ${quality} -T ${cpus} | pigz -p ${cpus} -c ${_ctx.input} > ${_ctx.output|ext}.map.gz"
  command += "gt.filter -i ${map} --max-levenshtein-error ${params.mismatches} -t ${cpus}| gt.filter --max-matches ${params.multimaps + 1} -t ${cpus} | pigz -p ${cpus} -c > ${_ctx.output|ext}.filter.map.gz"
  command += "pigz -p ${cpus} -dc ${filtered} | gem-2-sam -T ${cpus} -I ${index} -q ${quality} -l --expect-single-end-reads | awk '${awk_str}' | samtools view -@ ${cpus} -Sb - | samtools sort -@ ${cpus} - ${_ctx.output} > ${_ctx.output}"
  command += "samtools view -@ ${cpus} -bF256 ${bam}  > ${bam}_primary.bam"
  println "${command}"
}
return

/*process model {
  """
  run_spp.R -c=${bam} -rf -out=${prefix}.params.out -savp=${prefix}.pdf -p=${cpus}
  """
}

process peak {
  script:
  broad = { '--broad' if peak in broad_peaks else ''}
  command = ""
  command += "macs2 callpeak -t ${bam} -n ${prefix} --gsize hs -c ${chip_input} --nomodel --shiftsize=half_fragment_size ${broad}"
  """
}

process wiggle {

}*/
