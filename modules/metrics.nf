process NRF {
  input:
  tuple val(prefix), file(bam), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), stdout

  script:
  def cpus = task.cpus
  """
  samtools view -@${cpus} ${bam} | NRF.awk
  """
}

process FRiP {
  input:
  tuple val(prefix), file(bam), file(peak)

  output:
  tuple val(prefix), stdout

  script:
  def cpus = task.cpus
  """
  READS=\$(samtools view -F 260 -c ${bam})
  RiP=\$(bedtools intersect -abam ${bam} -b ${peak} -f 1.0 | samtools view - -c)
  awk -v rip=\$RiP -v reads=\$READS 'BEGIN{printf "%.4f", rip/reads}'
  """
}

workflow computeMetrics {
    take:
      bams
      peaks

    main:
      treatmentBams = bams
        .filter { it[3] != 'input' }
        .map {
            it[0..1] + it[3..-2]
        }
      FRiPbams = treatmentBams
        .mix(peaks.filter { it[2] != 'input' })
        .groupTuple(by: [0,2,3])
        .map {
          def bam = it[1].find { it.name =~ /bam$/ }
          def peaks = it[1].find { it.name =~ /narrowPeak$/ }
          [ it[0], bam, peaks ]
        }
      NRF(treatmentBams)
      FRiP(FRiPbams)
      metrics = NRF.out.cross(FRiP.out)
        .map { nrf, frip ->
          [nrf[0], nrf[1], frip[1]]
        }
    emit:
      metrics
}
