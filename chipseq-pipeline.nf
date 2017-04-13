params.dbFile = 'chipseq-pipeline.db'
params.genome = "${baseDir}/data/genome.fa"
params.genomeIndex = ''
params.genomeSize = 'hs'
params.fragmentLength = 200
params.help = false
params.index = "${baseDir}/data/index.tsv"
params.minMatchedBases = 0.8
params.mismatches = 2
params.multimaps = 10
params.qualityThreshold = 26
params.rescale = false
params.removeDuplicates = false
params.shift = false
params.zeroneConfidence  = 0.99
params.replicatePattern = '.[12]'

//print usage
if (params.help) {
    log.info ''
    log.info 'C H I P - N F ~ ChIP-seq Pipeline'
    log.info '---------------------------------'
    log.info 'Run ChIP-seq analyses on a set of data.'
    log.info ''
    log.info 'Usage: '
    log.info '    chipseq-pipeline.nf --index TSV_FILE --genome GENOME_FILE [OPTION]...'
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --index TSV_FILE                    Tab separted file containing information about the data.'
    log.info '    --genome GENOME_FILE                Reference genome file.'
    log.info '    --genome-index GENOME_INDEX_ FILE   Reference genome index file.'
    log.info '    --genome-size GENOME_SIZE           Reference genome size for MACS2 callpeaks. Must be one of' 
    log.info '                                        MACS2 precomputed sizes: hs, mm, dm, ce. (Default: hs)'
    log.info '    --mismatches MISMATCHES             Sets the maximum number/percentage of mismatches allowed for a read (Default: 2).'
    log.info '    --multimaps MULTIMAPS               Sets the maximum number of mappings allowed for a read (Default: 10).'
    log.info '    --min-matched-bases BASES           Sets the minimum number/percentage of bases that have to match with the reference (Default: 0.80).'
    log.info '    --quality-threshold THRESHOLD       Sets the sequence quality threshold for a base to be considered as low-quality (Default: 26).'
    log.info '    --fragment-length LENGTH            Sets the fragment length globally for all samples (Default: 200).'
    log.info '    --remove-duplicates                 Remove duplicate alignments instead of just flagging them (Default: false).'
    log.info '    --rescale                           Rescale peak scores to conform to the format supported by the'
    log.info '                                        UCSC genome browser (score must be <1000) (Default: false).'
    log.info '    --shift                             Move fragments ends and apply global extsize in peak calling. (Default: false).'
    log.info ''
    exit 1
}

pdb = file(params.dbFile)
pdb.write('')

////// Check input parameters //////
if (!params.genome) {
  exit 1, "Please specify a genome file"
}

if (!params.index) {
  exit 1, "Please specify the input table file"
}
////// End of input parameters check ////////

////// Print parameters ///////
log.info ''
log.info 'C H I P - N F ~ ChIP-seq Pipeline'
log.info '---------------------------------'
log.info ''
log.info "General parameters"
log.info '------------------'
log.info ''
log.info "Index File             : ${params.index}"
log.info "Genome File            : ${params.genome}"
log.info "Genome Index File      : ${params.genomeIndex ?: '-'}"
log.info "MACS2 Genome Size      : ${params.genomeSize}"
log.info "Global Fragment Length : ${params.fragmentLength}"
log.info "Database file          : ${params.dbFile}"
log.info "Remove Duplicates      : ${params.removeDuplicates}"
log.info "Shift                  : ${params.shift}"
log.info "Rescale Peaks          : ${params.rescale}"
log.info ''
log.info "Mapping parameters"
log.info '------------------'
log.info ''
log.info "Max Mismatches         : ${params.mismatches}"
log.info "Max Multimaps          : ${params.multimaps}"
log.info "Minimum Matched Bases  : ${params.minMatchedBases}"
log.info "Low Quality Threshold  : ${params.qualityThreshold}"
log.info ''
log.info "Zerone parameters"
log.info '-----------------'
log.info ''
log.info "Confidence Threshold   : ${params.zeroneConfidence}"
log.info ''

genome = file(params.genome)
index = file(params.index)

fastqs = Channel
.from(index.readLines())
.map { line ->
  def list = line.split()
  def mergeId = list[0]
  def id = list[1]
  def path = file(list[2])
  def controlId = list[3]
  def mark = list[4]
  def fragLen = list.size() == 6 ? list[5] as Integer : -1
  def message = '[INFO] '
  if (params.shift) {
      message += "Using global fragment length `${params.fragmentLength}` and compute shift size by "
  } else {
      message += "Using "
  }
  if ( fragLen != -1) {
    message += "fragment length `${fragLen}` for ${mergeId}"
  } else {
    message += "estimated fragment length for ${mergeId}"
  }
  log.info message
  def quality = fastq(path).qualityScore()
  [ mergeId, id, path, controlId, mark, fragLen, quality ]
}

if (!params.genomeIndex) {

  process index {

    input:
    file genome

    output:
    file "genome_index.gem" into GenomeIdx

    script:
    """
    sed 's/ .*//' ${genome} > genome_processed.fa
    gem-indexer -i genome_processed.fa \
                -o genome_index \
                -T ${task.cpus} \
                -m ${task.memory?.toBytes() ?: 'unlimited'} \
    && rm genome_processed.fa
    """
  }

} else {
    GenomeIdx = Channel.fromPath(params.genomeIndex)
}

process fastaIndex {
  input:
  file genome

  output:
  file "${genome}.fai" into chromSizesNarrowPeakCall, chromSizesNarrowPeakCallNoInput, chromSizesBroadPeakCall, chromSizesBroadPeakCallNoInput, chromSizesPileupSignalTracks, chromSizesFeSignalTracks, chromSizesPvalSignalTracks

  script:
  """
  samtools faidx ${genome}
  """
}

process mapping {
  input:
  file index from GenomeIdx.val
  set mergeId, prefix, file(fastq), controlId, mark, fragLen, quality from fastqs

  output:
  set mergeId, prefix, file("${prefix}_primary.bam"), controlId, mark, fragLen, val('Alignments') into bams

  script:
  def cpus = task.cpus
  def memory = task.memory
  def readGroup = "ID=${prefix},SM=${mergeId}"
  def cat = fastq.name.endsWith('.gz') ? 'zcat' : 'cat'
  def awk_str = 'BEGIN{OFS=FS="\\t"}$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and($2,a[i])>0){$2=xor($2,a[i])}}}{print}'
  """
  ${cat} ${fastq} \
  | gem-mapper -I ${index} \
               -m ${params.mismatches} \
               -e ${params.mismatches} \
               --min-matched-bases ${params.minMatchedBases} \
               -d ${params.multimaps} \
               -D 0 \
               --gem-quality-threshold ${params.qualityThreshold} \
               -q offset-${quality} \
               -T ${cpus} \
  | gem-2-sam -T ${cpus} \
              -I ${index} \
              -q offset-${quality} \
              -l \
              --expect-single-end-reads \
              --read-group ${readGroup} \
  | awk '${awk_str}' \
  | samtools view -@ ${cpus} \
                  -SbF256 \
                  - \
  | samtools sort -@ ${cpus} \
                  - \
                  ${prefix}_primary
  """
}
// Merge or rename bam
singleBams = Channel.create()
groupedBams = Channel.create()

bams.groupTuple(by: [0,3,4])
.choice(singleBams, groupedBams) {
  it[2].size() > 1 ? 1 : 0
}

process mergeBam {

    input:
    set mergeId, prefix, file(bam), controlId, mark, fragLen, view from groupedBams

    output:
    set mergeId, prefix, file("${mergeId}_primary.bam"), controlId, mark, fragLen, view into mergedBams

    script:
    def cpus = task.cpus
    def prefix = prefix.sort().join(':')
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


singleBams
.mix(mergedBams)
.map { mergeId, prefix, bam, controlId, mark, fragLen, view ->
  [ mergeId, bam, controlId, mark, fragLen, view].flatten()
}
.into { bamsMarkDup }


process markDup {
  
  input:
  set prefix, file(bam), controlId, mark, fragLen, view from bamsMarkDup

  output:
  set prefix, file("${bam.baseName}_picard.bam"), controlId, mark, fragLen, view into bamsMarked, originalBams, bamsReadCount, modelBams
  set prefix, file("${prefix}.picard.metrics"), controlId, mark, fragLen, val("picardMetrics") into picardStats

  script:
  def sorted = true
  def rmDup = params.removeDuplicates
  def mem = "${task.memory?.toMega() ?: 1024}m"
  """
  java -Xmx${mem} -jar `which picard.jar` MarkDuplicates INPUT=${bam} \
                                               OUTPUT=${bam.baseName}_picard.bam \
                                               METRICS_FILE=${prefix}.picard.metrics \
                                               VALIDATION_STRINGENCY=LENIENT \
                                               ASSUME_SORTED=${sorted} \
                                               REMOVE_DUPLICATES=${rmDup}
  """
}

process readCount {

//  when:
//  controlId != '-'

  input:
  set prefix, file(bam), controlId, mark, fragLen, view from bamsReadCount

  output:
  set prefix, stdout into bamsReads

  script:
  """
  samtools view -F 256 -c ${bam} | tr -d '\n'
  """

}

process model {
  when:
  fragLen == -1

  input:
  set prefix, file(bam), controlId, mark, fragLen, view from modelBams

  output:
  set prefix, file("${prefix}.params.out") into modelParams

  script:
  def cpus = task.cpus
  """
  Rscript \$(which run_spp.R) -c=${bam} \
                              -rf \
                              -out=${prefix}.params.out \
                              -savp=${prefix}.pdf \
                              -p=${cpus}
  """
}

originalBams.cross(bamsReads)
.map { bams, reads ->
  bams + reads[-1..-1]
}.mix(modelParams)
.groupBy()
.flatMap()
.map { it ->
  switch (it.value.size()) {
    case 1:
      it.value[0]
      break
    case 2:
      def bams = it.value.find { it[1] =~ /bam/ }
      def paramFile = it.value.find { it[1] =~ /params/ }[1]
      def fragLen = paramFile.text.split()[2].split(',')[0] as Integer
      bams[0..-4] + [fragLen] + bams[-2..-1]
      break
  }
}.tap{ allBams }
.map { prefix, bam, controlId, mark, fragLen, view, readCount->
  [prefix, bam, controlId, mark, fragLen, view]
}.tap { bamResults }
.filter {
  it[3] != 'input'
}.into { bams4NRF; bams4FRiP }

// separate bams and inputs
treatBams = Channel.create()
controlBams = Channel.create()
allBams.choice(treatBams, controlBams) {
    it[3] == 'input' ? 1 : 0
}

// get bams with no control
treatBams.tap { inputBams }
.filter {
  it[2] == '-'
}
.map { prefix, bam, controlId, mark, fragLen, view, readCount ->
  [ prefix, bam, mark, fragLen, view ]
}
.into { bamsNarrowPeakCallNoInput; bamsBroadPeakCallNoInput }

// cross bams and controls
controlBams
.tap { controlBams4Signal }
.filter {
  it[2] != '-'
}
.phase(inputBams) { it[2] }
.tap { crossedBams } 
.map { c, t ->
  [t[0], t[1], c[1], t[3], t[4], t[5]]
}
.tap { bamsNarrowPeakCall; bamsBroadPeakCall }
.map { replicateId, bam, control, mark, fragLen, view ->
  sampleId = replicateId.replaceAll(/${params.replicatePattern}$/,'')
  [sampleId, bam, control, mark, view]
}
.groupTuple(by:[0,3,4])
.into {bamsZerone}

process narrowPeakCall {
  
  input:
  file chromSizes from chromSizesNarrowPeakCall.val
  set prefix, file(bam), file(control), mark, fragLen, view from bamsNarrowPeakCall

  output:
  set prefix, file("peakOut/${prefix}_peaks.narrowPeak"), mark, fragLen, val("narrowPeak") into narrowPeakFiles, narrowPeakFiles4FRiP
  set prefix, file("peakOut/${prefix}*.bdg"), mark, fragLen, val("pileupBedGraphs") into pileupBedGraphFiles, pileupBedGraphFilesPileupSignalTracks, pileupBedGraphFilesFeSignalTracks
  
  script:
  def extSize = params.shift ? params.fragmentLength : fragLen
  def shiftSize = params.shift ? Math.round((fragLen - params.fragmentLength) / 2) : 0
  """
  # narrow peaks and preliminary signal tracks
  macs2 callpeak -t ${bam} -c ${control} -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 \
                 --nomodel --shift=${shiftSize} --extsize=${extSize} \
                 --keep-dup all -B --SPMR
  """
}

process narrowPeakCallNoInput {
  
  input:
  file chromSizes from chromSizesNarrowPeakCallNoInput.val
  set prefix, file(bam), mark, fragLen, view from bamsNarrowPeakCallNoInput.mix(controlBams4Signal.map { [ it[0], it[1], it[3], it[4], it[5] ] })

  output:
  set prefix, file("peakOut/${prefix}_peaks.narrowPeak"), mark, fragLen, val("narrowPeak") into narrowPeakFilesNoInput, narrowPeakFiles4FRiPNoInput
  set prefix, file("peakOut/${prefix}*.bdg"), mark, fragLen, val("pileupBedGraphs") into pileupBedGraphFilesNoInput, pileupBedGraphFilesPileupSignalTracksNoInput, pileupBedGraphFilesFeSignalTracksNoInput
  
  script:
  def extSize = params.shift ? params.fragmentLength : fragLen
  def shiftSize = params.shift ? Math.round((fragLen - params.fragmentLength) / 2) : 0
  """
  # narrow peaks and preliminary signal tracks
  macs2 callpeak -t ${bam} -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 \
                 --nomodel --shift=${shiftSize} --extsize=${extSize} \
                 --keep-dup all -B --SPMR
  """
}

crossedBams.map{ c, t ->
  [t[0], t[4], c[4]]
}.cross(pileupBedGraphFiles)
.map { r, s ->
  def (treat, control) = r[1..-1] as long[]
  def count = treat < control ? treat : control
  [s[0], s[1], s[2], s[3], count/1000000, s[4]]
}
.into{ pileupBedGraphFilesPvalSignalTracks }

process broadPeakCall {
  
  input:
  file chromSizes from chromSizesBroadPeakCall.val
  set prefix, file(bam), file(control), mark, fragLen, view from bamsBroadPeakCall

  output:
  set prefix, file("peakOut/${prefix}_peaks.broadPeak"), mark, fragLen, val("broadPeak") into broadPeakFiles
  set prefix, file("peakOut/${prefix}_peaks.gappedPeak"), mark, fragLen, val("gappedPeak") into gappedPeakFiles

  script:
  def extSize = params.shift ? params.fragmentLength : fragLen
  def shiftSize = params.shift ? Math.round((fragLen - params.fragmentLength) / 2) : 0
  """
  # Broad and Gapped Peaks
  macs2 callpeak -t ${bam} -c ${control} -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 --broad \
                 --nomodel --shift=${shiftSize} --extsize=${extSize} \
                 --keep-dup all
  """
}

process broadPeakCallNoInput {
  
  input:
  file chromSizes from chromSizesBroadPeakCallNoInput.val
  set prefix, file(bam), mark, fragLen, view from bamsBroadPeakCallNoInput

  output:
  set prefix, file("peakOut/${prefix}_peaks.broadPeak"), mark, fragLen, val("broadPeak") into broadPeakFilesNoInput
  set prefix, file("peakOut/${prefix}_peaks.gappedPeak"), mark, fragLen, val("gappedPeak") into gappedPeakFilesNoInput

  script:
  def extSize = params.shift ? params.fragmentLength : fragLen
  def shiftSize = params.shift ? Math.round((fragLen - params.fragmentLength) / 2) : 0
  """
  # Broad and Gapped Peaks
  macs2 callpeak -t ${bam} -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 --broad \
                 --nomodel --shift=${shiftSize} --extsize=${extSize} \
                 --keep-dup all
  """
}

// Collect peak files
narrowPeakFiles
  .mix( broadPeakFiles )
  .mix( gappedPeakFiles )
  .mix( narrowPeakFilesNoInput )
  .mix( broadPeakFilesNoInput )
  .mix( gappedPeakFilesNoInput )
.into{ allPeakFiles; peakCallResults; peakCallResults4FRiP }

process rescalePeaks {
 
  when:
  params.rescale

  input:
  set prefix, file(peaks), mark, fragLen, view from allPeakFiles

  output:
  set prefix, file("${peaks.name}.rescaled"), mark, fragLen, val("rescaled${view.capitalize()}") into rescaledPeakFiles


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
  file chromSizes from chromSizesPileupSignalTracks.val
  set prefix, file(bedGraphs), mark, fragLen, view from pileupBedGraphFilesPileupSignalTracks.mix(pileupBedGraphFilesPileupSignalTracksNoInput)

  output:
  set prefix, file("${prefix}.pileup_signal.bw"), mark, fragLen, val("pileupSignal") into pileupSignalFiles, pileupSignalFilesResults
  
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
  file chromSizes from chromSizesFeSignalTracks.val
  set prefix, file(bedGraphs), mark, fragLen, view from pileupBedGraphFilesFeSignalTracks.mix(pileupBedGraphFilesFeSignalTracksNoInput)

  output:
  set prefix, file("${prefix}.fc_signal.bw"), mark, fragLen, val("fcSignal") into feSignalFiles

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
  file chromSizes from chromSizesPvalSignalTracks.val
  set prefix, file(bedGraphs), mark, fragLen, sval, view from pileupBedGraphFilesPvalSignalTracks
 
  output:
  set prefix, file("${prefix}.pval_signal.bw"), mark, fragLen, val("pvalueSignal") into pvalSignalFiles
  
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

process zerone {

  input:
  set prefix, file(bam), file(control), mark, view from bamsZerone

  output:
  set prefix, file("${prefix}_zerone.01"), mark, val("zeroneMatrix") into zeroneMatrixFiles
  set prefix, file("${prefix}_zerone.bed"), mark, val("zeroneBed") into zeroneBedFiles
  set prefix, file("${prefix}_zerone_merged.bed"), mark, val("zeroneMergedBed") into zeroneMergedBedFiles

  script:
  def awkScaleMergedBed = '$0~/^#/ || $NF=$NF*1000'
  def awkMatrix2Bed = '$0~/^#/ || $0=$1 OFS $2 OFS $3 OFS $NF*1000'
  """
  zerone -c ${params.zeroneConfidence} -0 ${control.join(",")} -1 ${bam.join(",")} > ${prefix}_zerone.01
  awk -F"\\t" '${awkMatrix2Bed}' OFS="\\t" ${prefix}_zerone.01 > ${prefix}_zerone.bed
  zerone -c ${params.zeroneConfidence} -l -0 ${control.join(",")} -1 ${bam.join(",")} | awk -F"\\t" '${awkScaleMergedBed}' OFS="\\t" > ${prefix}_zerone_merged.bed
  """
}

process NRF {
  input:
  set prefix, file(bam), controlId, mark, view from bams4NRF

  output:
  set prefix, stdout into NRFBams

  script:
  def cpus = task.cpus
  """
  samtools view -@${cpus} ${bam} | NRF.awk 
  """
}

input4FRiP = bams4FRiP.map { prefix, bam, controlId, mark, fragLen, view ->
    [prefix, bam, mark, fragLen, view]
}.mix(narrowPeakFiles4FRiP)
.mix(narrowPeakFiles4FRiPNoInput)
.filter { it[2] != 'input' }
.groupTuple(by: [0,2,3])
.map { prefix, files, controlId, mark, views ->
  def bam = files.find { it.name =~ /bam$/ }
  def peaks = files.find { it.name =~ /narrowPeak$/ }
  [ prefix, bam, peaks ]
}

process FRiP {

  input:
  set prefix, file(bam), file(peak) from input4FRiP

  output:
  set prefix, stdout into FRiPBams

  script:
  def cpus = task.cpus
  """
  READS=\$(samtools view -F 260 -c ${bam})
  RiP=\$(bedtools intersect -abam ${bam} -b ${peak} -f 1.0 | samtools view - -c)
  echo "\$RiP/\$READS" | bc -l | awk '{printf "%.4f", \$0}'
  """
}

metrics = NRFBams.cross(FRiPBams)
.map { nrf, frip ->
 [nrf[0], nrf[1], frip[1]] 
}

metrics.cross(
  bamResults
  .map { prefix, bam, control, mark, fragLen, view ->
    [ prefix, bam, mark, fragLen, view ]
  }
  .tap { resultBams }
  .mix(peakCallResults, pileupSignalFiles, feSignalFiles, pvalSignalFiles)
).map { qc, result ->
    result + qc[1..-1]
}.mix(
  resultBams.mix(pileupSignalFilesResults)
  .filter { it[2] == 'input' }
  .map {
    it + [ '-', '-' ]
  }
).mix(
  zeroneMatrixFiles.mix(zeroneBedFiles).mix(zeroneMergedBedFiles)
  .map { prefix, path, mark, view ->
    [prefix, path, mark, '-', view, '-', '-']
  }
)
.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) { prefix, path, mark, fragLen, view, nrf, frip ->
    [ prefix, path, mark, fragLen, view, nrf, frip ].join("\t")
}


workflow.onComplete {
    log.info ""
    log.info "-----------------------"
    log.info "Pipeline run completed."
    log.info "-----------------------"
}
