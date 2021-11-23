include { fastaIndex } from './ref'

process narrowPeak {
  input:
  file chromSizes
  val globalFragmentLength
  tuple val(prefix), file(bam), file(control), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("peakOut/${prefix}_peaks.narrowPeak"), val(mark), val(fragLen), val("narrowPeak"), emit: narrowPeaks
  tuple val(prefix), file("peakOut/${prefix}*.bdg"), val(mark), val(fragLen), val("pileupBedGraphs"), emit: pileupBedGraphs

  script:
  def extSize = params.shift ? globalFragmentLength : fragLen
  def shiftSize = params.shift ? Math.round((fragLen - globalFragmentLength) / 2) : 0
  """
  # narrow peaks and preliminary signal tracks
  macs2 callpeak -t ${bam} -c ${control} -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 \
                 --nomodel --shift=${shiftSize} --extsize=${extSize} \
                 --keep-dup all -B --SPMR --tempdir ${params.tmpDirMACS2}
  """
}

process narrowPeakNoInput {
  input:
  file chromSizes
  val globalFragmentLength
  tuple val(prefix), file(bam), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("peakOut/${prefix}_peaks.narrowPeak"), val(mark), val(fragLen), val("narrowPeak"), emit: narrowPeaks
  tuple val(prefix), file("peakOut/${prefix}*.bdg"), val(mark), val(fragLen), val("pileupBedGraphs"), emit: pileupBedGraphs

  script:
  def extSize = params.shift ? globalFragmentLength : fragLen
  def shiftSize = params.shift ? Math.round((fragLen - globalFragmentLength) / 2) : 0
  """
  # narrow peaks and preliminary signal tracks
  macs2 callpeak -t ${bam} -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 \
                 --nomodel --shift=${shiftSize} --extsize=${extSize} \
                 --keep-dup all -B --SPMR --tempdir ${params.tmpDirMACS2}
  """
}

process broadPeak {
  input:
  file chromSizes
  val globalFragmentLength
  tuple val(prefix), file(bam), file(control), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("peakOut/${prefix}_peaks.broadPeak"), val(mark), val(fragLen), val("broadPeak"), emit: broadPeaks
  tuple val(prefix), file("peakOut/${prefix}_peaks.gappedPeak"), val(mark), val(fragLen), val("gappedPeak"), emit: gappedPeaks

  script:
  def extSize = params.shift ? globalFragmentLength : fragLen
  def shiftSize = params.shift ? Math.round((fragLen - globalFragmentLength) / 2) : 0
  """
  # Broad and Gapped Peaks
  macs2 callpeak -t ${bam} -c ${control} -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 --broad \
                 --nomodel --shift=${shiftSize} --extsize=${extSize} \
                 --keep-dup all --tempdir ${params.tmpDirMACS2}
  """
}

process broadPeakNoInput {
  input:
  file chromSizes
  val globalFragmentLength
  tuple val(prefix), file(bam), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("peakOut/${prefix}_peaks.broadPeak"), val(mark), val(fragLen), val("broadPeak"), emit: broadPeaks
  tuple val(prefix), file("peakOut/${prefix}_peaks.gappedPeak"), val(mark), val(fragLen), val("gappedPeak"), emit: gappedPeaks

  script:
  def extSize = params.shift ? globalFragmentLength : fragLen
  def shiftSize = params.shift ? Math.round((fragLen - globalFragmentLength) / 2) : 0
  """
  # Broad and Gapped Peaks
  macs2 callpeak -t ${bam} -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 --broad \
                 --nomodel --shift=${shiftSize} --extsize=${extSize} \
                 --keep-dup all --tempdir ${params.tmpDirMACS2}
  """
}

process rescalePeaks {
  input:
  tuple val(prefix), file(peaks), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("${peaks.name}.rescaled"), val(mark), val(fragLen), val("rescaled${view.capitalize()}")


  script:
  def rescale_awk_str = 'BEGIN{FS=OFS="\\t";min=1e20;max=-1}'
  rescale_awk_str += 'NR==FNR&&NF!=0{min>$5?min=$5:min=min;max<$5?max=$5:max=max;next}'
  rescale_awk_str += 'NF!=0{n=$5;x=10;y=1000;$5=int(((n-min)*(y-x)/(max-min))+x);print}'
  """
  # rescale peaks on 10-1000 scale
  awk '${rescale_awk_str}' ${peaks} > ${peaks.name}.rescaled
  """
}

process pileupSignalTracks {
  input:
  file chromSizes
  tuple val(prefix), file(bedGraphs), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("${prefix}.pileup_signal.bw"), val(mark), val(fragLen), val("pileupSignal")

  script:
  def treat = bedGraphs instanceof nextflow.util.BlankSeparatedList ? bedGraphs.find { it =~ /treat/ } : bedGraphs
  """
  # pileup signal tracks
  slopBed -i ${treat} -g ${chromSizes} -b 0 \
  | bedClip stdin ${chromSizes} ${prefix}.pileup_signal.bedgraph
  bedGraphToBigWig ${prefix}.pileup_signal.bedgraph ${chromSizes} ${prefix}.pileup_signal.bw
  """
}

process feSignalTracks {
  when:
  bedGraphs instanceof nextflow.util.BlankSeparatedList && mark != 'input'

  input:
  file chromSizes
  tuple val(prefix), file(bedGraphs), val(mark), val(fragLen), val(view)

  output:
  tuple val(prefix), file("${prefix}.fc_signal.bw"), val(mark), val(fragLen), val("fcSignal")

  script:
  def treat = bedGraphs.find { it =~ /treat/ }
  def control = bedGraphs.find { it =~ /control/ }
  """
  # Fold enrichment signal tracks
  macs2 bdgcmp -t ${treat} \
               -c ${control} --outdir . \
               -o ${prefix}_FE.bdg -m FE
  slopBed -i ${prefix}_FE.bdg -g ${chromSizes} -b 0 \
  | bedClip stdin ${chromSizes} ${prefix}.fc.signal.bedgraph
  bedGraphToBigWig ${prefix}.fc.signal.bedgraph ${chromSizes} ${prefix}.fc_signal.bw
  """
}

process pvalSignalTracks {
  when:
  bedGraphs instanceof nextflow.util.BlankSeparatedList && mark != 'input'

  input:
  file chromSizes
  tuple val(prefix), file(bedGraphs), val(mark), val(fragLen), val(sval), val(view)

  output:
  tuple val(prefix), file("${prefix}.pval_signal.bw"), val(mark), val(fragLen), val("pvalueSignal")

  script:
  def treat = bedGraphs.find { it =~ /treat/ }
  def control = bedGraphs.find { it =~ /control/ }
  """
  # -log10(p-value) signal tracks
  macs2 bdgcmp -t ${treat} \
               -c ${control} --outdir . \
               -o ${prefix}_ppois.bdg -m ppois -S ${sval}
  slopBed -i ${prefix}_ppois.bdg -g ${chromSizes} -b 0 \
  | bedClip stdin ${chromSizes} ${prefix}.pval.signal.bedgraph
  bedGraphToBigWig ${prefix}.pval.signal.bedgraph ${chromSizes} ${prefix}.pval_signal.bw
  """
}

workflow callPeaks {
    take:
      genome
      globalFragmentLength
      bams
    main:
      controlBams = bams.filter { it[3] == 'input' }
      treatmentBams = bams.filter { it[2] != '-' && it[3] != 'input' }
      crossedBams = controlBams.cross( treatmentBams ) { it[2] }
      treatmentBamsWithInput = crossedBams
        .map { control, treatment ->
            // Add control BAM file
            treatment[0..1] + control[1..1] + treatment[3..-2]
        }
      treatmentBamsNoInput = bams
        .filter { it[2] == '-' }
        .map {
            // Remove 'input' column
            it[0..1] + it[3..-2]
        }
      treatmentBamsNoInputAndControls = controlBams
        .map {
            // Remove 'input' column
            it[0..1] + it[3..-2]
        }.mix(treatmentBamsNoInput)
      fastaIndex(genome)
      narrowPeak(fastaIndex.out, globalFragmentLength, treatmentBamsWithInput)
      narrowPeakNoInput(fastaIndex.out, globalFragmentLength, treatmentBamsNoInputAndControls)
      broadPeak(fastaIndex.out, globalFragmentLength, treatmentBamsWithInput)
      broadPeakNoInput(fastaIndex.out, globalFragmentLength, treatmentBamsNoInput)

      allPeaks = narrowPeak.out.narrowPeaks
        .mix(narrowPeakNoInput.out.narrowPeaks)
        .mix(broadPeak.out.broadPeaks)
        .mix(broadPeak.out.gappedPeaks)
        .mix(broadPeakNoInput.out.broadPeaks)
        .mix(broadPeakNoInput.out.gappedPeaks)
      if (params.rescale)
        rescalePeaks(allPeaks)

      allBedGraphs = narrowPeak.out.pileupBedGraphs
        .mix(narrowPeakNoInput.out.pileupBedGraphs)

      pileupSignalTracks(fastaIndex.out, allBedGraphs)
      feSignalTracks(fastaIndex.out, allBedGraphs)

      allBedGraphsScaleFactor = crossedBams.map { c, t ->
        [t[0], t[6], c[6]]
      }.cross(allBedGraphs)
      .map { r, s ->
        def (treat, control) = r[1..-1] as long[]
        def count = treat < control ? treat : control
        [s[0], s[1], s[2], s[3], count/1000000, s[4]]
      }
      pvalSignalTracks(fastaIndex.out, allBedGraphsScaleFactor)
      allBigWigs = pileupSignalTracks.out
        .mix(feSignalTracks.out)
        .mix(pvalSignalTracks.out)

    emit:
      params.rescale ? rescalePeaks.out : allPeaks
      allBigWigs
}
