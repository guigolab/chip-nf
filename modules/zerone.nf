process zerone {
  input:
  tuple val(prefix), file(bam), file(control), val(mark), val(view)

  output:
  tuple val(prefix), file("${prefix}_zerone.01"), val(mark), val("zeroneMatrix")
  tuple val(prefix), file("${prefix}_zerone.bed"), val(mark), val("zeroneBed")
  tuple val(prefix), file("${prefix}_zerone_merged.bed"), val(mark), val("zeroneMergedBed")

  script:
  def awkScaleMergedBed = '$0~/^#/ || $NF=$NF*1000'
  def awkMatrix2Bed = '$0~/^#/ || $0=$1 OFS $2 OFS $3 OFS $NF*1000'
  """
  zerone -c ${params.zeroneMinConfidence} -0 ${control.join(",")} -1 ${bam.join(",")} > ${prefix}_zerone.01
  awk -F"\\t" '${awkMatrix2Bed}' OFS="\\t" ${prefix}_zerone.01 > ${prefix}_zerone.bed
  zerone -c ${params.zeroneMinConfidence} -l -0 ${control.join(",")} -1 ${bam.join(",")} | awk -F"\\t" '${awkScaleMergedBed}' OFS="\\t" > ${prefix}_zerone_merged.bed
  """
}

workflow callPeaks {
  take: bams
  main:
    controlBams = bams.filter { it[3] == 'input' }
    treatmentBams = bams.filter { it[2] != '-' && it[3] != 'input' }
    crossedBams = controlBams.cross( treatmentBams ) { it[2] }
    zeroneBams = crossedBams.map { c, t ->
      // Join replicates
      sampleId = t[0].replaceAll(/${params.replicatePattern}$/,'')
      [ sampleId, t[1], c[1], t[3], t[5] ]
    }
    .groupTuple(by:[0,3,4], sort: {it.baseName})
    .map{
      it[0..1] + [it[2].unique()] + it[3..-1]
    }
    (matrix, bed, mergedBed) = zerone(zeroneBams)
  emit:
    matrix
    bed
    mergedBed
}
