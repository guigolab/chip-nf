process merge {
    input:
    tuple val(mergeId), val(prefix), file(bam), val(controlId), val(mark), val(fragLen), val(view)

    output:
    tuple val(mergeId), val(prefix), file("${mergeId}_primary.bam"), val(controlId), val(mark), val(fragLen), val(view)

    script:
    def cpus = task.cpus
    prefix = prefix.sort().join(':')
    """
    (
      samtools view -H ${bam} | grep -v '@RG';
      for f in ${bam}; do
        samtools view -H \$f | grep '@RG';
      done
    ) > header.txt && \
    samtools merge -@ ${cpus} \
                   -h header.txt \
                   ${mergeId}_primary.bam \
                   ${bam}
    """
}

process markDup {
  input:
  tuple val(prefix), file(bam), val(controlId), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("${bam.baseName}_picard.bam"), val(controlId), val(mark), val(fragLen), val(view), emit: bams
  tuple val(prefix), file("${prefix}.picard.metrics"), val(controlId), val(mark), val(fragLen), val("picardMetrics"), emit: metrics

  script:
  def sorted = true
  def rmDup = params.removeDuplicates
  def mem = "${task.memory?.toMega() ?: 1024}m"
  """
  java -Xmx${mem} -jar `which picard.jar` MarkDuplicates \
                            INPUT=${bam} \
                            OUTPUT=${bam.baseName}_picard.bam \
                            METRICS_FILE=${prefix}.picard.metrics \
                            VALIDATION_STRINGENCY=LENIENT \
                            ASSUME_SORTED=${sorted} \
                            REMOVE_DUPLICATES=${rmDup}
  """
}

process readCount {
  input:
  tuple val(prefix), file(bam), val(controlId), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), stdout

  script:
  """
  samtools view -F 256 -c ${bam} | tr -d '\\n'
  """

}

process model {
  input:
  tuple val(prefix), file(bam), val(controlId), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("${prefix}.params.out")

  script:
  def cpus = task.cpus
  """
  Rscript \$(which run_spp.R) -c=${bam} \
                              -rf \
                              -out=${prefix}.params.out \
                              -savp=${prefix}.pdf \
                              ${ cpus > 1 ? "-p=${cpus}" : '' }
  """
}

workflow mergeBams {
    take: bams
    main:
      // group by 'mergeId','controlId','mark','fragLen','view'
      groupedBams = bams.groupTuple(by: [0,3,4,5,6])

      bams2merge = groupedBams.filter { it[2].size() > 1 }
      mergedBams = merge(bams2merge)

      allBams = groupedBams
        .filter { it[2].size() == 1 }
        .mix(mergedBams)
        .map {
          // Remove 'prefix'
          it[0..0] + it[2..-1]
        }
    emit: allBams
}

workflow inferFragLen {
    take: bams
    main:
      model(bams.filter { !it[4] })
      allBams = bams.mix(model.out)
        .groupTuple()
        .map {
            if (it[1].size() == 1)
                return it.flatten()
            def bams = it[1].find { it =~ /bam/ }
            def paramFile = it[1].find { it =~ /params/ }
            def fragLen = paramFile.text.split()[2].split(',')[0] as Integer
            it[1] = bams; it[4] = fragLen; it.flatten()
        }
    emit:
      allBams
}

workflow processBams {
  take: bams
  main:
    markDup(bams)
    readCount(bams)
    inferFragLen(bams)
    allBams = inferFragLen.out
      .cross(markDup.out.bams)
      .map { bams, marked ->
        [ bams[0], marked[1] ] + bams[2..-1]
      }
      .cross(readCount.out)
      .map { bams, reads ->
        bams + reads[-1..-1]
      }
    emit:
      allBams
}
