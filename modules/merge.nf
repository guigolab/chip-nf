process mergeBam {

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

workflow mergeBams {
    take: bams
    main:
      // group by 'mergeId','controlId','mark','fragLen','view'
      groupedBams = bams.groupTuple(by: [0,3,4,5,6])

      bams2merge = groupedBams.filter { it[2].size() > 1 }
      mergedBams = mergeBam(bams2merge)

      allBams = groupedBams
        .filter { it[2].size() == 1 }
        .mix(mergedBams)
        .map {
          // Remove 'prefix'
          it[0..0] + it[2..-1]
        }
    emit: allBams
}

