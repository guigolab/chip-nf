params.dbFile = 'chipseq-pipeline.db'
params.genome = ''
params.genomeIndex = ''
params.genomeSize = 'hs'
params.help = false
params.index = ''
params.mismatches = 2
params.multimaps = 10
params.rescale = false

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
    log.info '    --mismatches N_MISMATCHES           Allow max N_MISMATCHES error events for a read (Default: 2).'
    log.info '    --multimaps N_MULTIMAPS             Allow max N_MULTIMAPS mappings for a read (Default: 10).'
    log.info '    --rescale                           Rescale peak scores to conform to the format supported by the '
    log.info '                                        UCSC genome browser (score must be <1000) (Default: false).'
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

genome = file(params.genome)
index = file(params.index)

fastqs = Channel
.from(index.readLines())
.map { line ->
  list = line.split()
  mergeId = list[0]
  id = list[1]
  path = file(list[2])
  controlId = list[3]
  mark = list[4]
  quality = fastq(path).qualityScore()
  [ mergeId, id, path, controlId, mark, quality ]
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
  file "${genome}.fai" into chromSizes

  script:
  """
  samtools faidx ${genome}
  """
}

process mapping {
  input:
  file index from GenomeIdx.val
  set mergeId, prefix, file(fastq), controlId, mark, quality from fastqs

  output:
  set mergeId, prefix, file("${prefix}_primary.bam"), controlId, mark, view into bams

  script:
  cpus = task.cpus
  memory = task.memory
  readGroup = "ID=${prefix},SM=${mergeId}"
  view = "Alignments"
  cat = fastq.name.endsWith('.gz') ? 'zcat' : 'cat'
  awk_str = 'BEGIN{OFS=FS="\\t"}$0!~/^@/{split(\"1_2_8_32_64_128\",a,\"_\");for(i in a){if(and($2,a[i])>0){$2=xor($2,a[i])}}}{print}'
  """
  ${cat} ${fastq} \
  | gem-mapper -I ${index} \
               -q offset-${quality} \
               -T ${cpus} \
  | pigz -p ${cpus} \
         -c \
  > ${prefix}.map.gz
  gt.filter -i ${prefix}.map.gz \
            --max-levenshtein-error ${params.mismatches} \
            -t ${cpus}\
  | gt.filter --max-maps ${params.multimaps} \
              -t ${cpus} \
  | pigz -p ${cpus} 
         -c \
  > ${prefix}.filter.map.gz
  pigz -p ${cpus} \
       -dc ${prefix}.filter.map.gz \
  | gem-2-sam -T ${cpus} \
              -I ${index} \
              -q offset-${quality} \
              -l \
              --expect-single-end-reads \
              --read-group ${readGroup} \
  | awk '${awk_str}' \
  | samtools view -@ ${cpus} \
                  -Sb \
                  - \
  | samtools sort -@ ${cpus} \
                  - \
                  ${prefix}
  
  samtools view -@ ${cpus} \
                -bF256 \
                ${prefix}.bam \
  > ${prefix}_primary.bam
  """
}

// Merge or rename bam
singleBam = Channel.create()
groupedBam = Channel.create()

bams.groupTuple(by: [0,3,4])
.choice(singleBam, groupedBam) {
  it[2].size() > 1 ? 1 : 0
}

process mergeBam {

    input:
    set mergeId, prefix, file(bam), controlId, mark, view from groupedBam

    output:
    set mergeId, prefix, file("${mergeId}.bam"), controlId, mark, view into mergedBam

    script:
    cpus = task.cpus
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


bams = singleBam
.mix(mergedBam)
.map { mergeId, prefix, bam, controlId, mark, view ->
  [ mergeId, bam, controlId, mark, view].flatten()
}

// separate bams and inputs
treat = Channel.create()
control = Channel.create()
bams.choice(treat, control) {
    it[3] == 'input' ? 1 : 0
}

process model {
  input:
  set prefix, file(bam), controlId, mark, view from treat

  output:
  set prefix, file("${prefix}.params.out") into modelParams
  set prefix, file(bam), controlId, file("${prefix}.params.out"), mark, view into modelBams

  script:
  cpus = task.cpus
  """
  Rscript \$(which run_spp.R) -c=${bam} \
                              -rf \
                              -out=${prefix}.params.out \
                              -savp=${prefix}.pdf \
                              -p=${cpus}
  """
}

(bams, results, bams4NRF, bams4FRiP ) = modelBams.map { prefix, bam, controlId, paramFile, mark, view ->
  fragLen = paramFile.text.split()[2].split(',')[0]
  [prefix, bam, controlId, mark, fragLen, view]
}.into(4)

// get bams with no control
bams.tap { allBams }
.filter {
  it[2] == '-'
}.map {
  [it[0], it[1], it[3], it[4], it[5]]
}.tap { bamsNoInput }

// cross bams and controls
bamsWithInput = control.filter {
  it[2] != '-'
}
.cross(allBams) { it[2] }.map { c, t ->
  [t[0], t[1], c[1], t[3], t[4], t[5]]
}

(chromSizesForInput, chromSizesForNoInput) = chromSizes.into(2)

process peakCall {
  input:
  file chromSizes from chromSizesForInput.val
  set prefix, file(bam), file(control), mark, fragLen, view from bamsWithInput

  output:
  set prefix, file("${prefix}*{Peak,bw}"), mark, fragLen into peakCallWithInputResults
  // set prefix, file("peakOut/${prefix}_peaks.narrowPeak"), mark, fragLen, val("narrowPeak") into peakCallWithInputResults
  // set prefix, file("peakOut/${prefix}_peaks.broadPeak"), mark, fragLen, val("broadPeak") into peakCallWithInputResults
  // set prefix, file("peakOut/${prefix}_peaks.gappedPeak"), mark, fragLen, val("gappedPeak") into peakCallWithInputResults
  // set prefix, file("peakOut/${prefix}.pileup_signal.bw"), mark, fragLen, val("pileupSignal") into peakCallWithInputResults
  // set prefix, file("peakOut/${prefix}.fc_signal.bw"), mark, fragLen, val("fcSignal") into peakCallWithInputResults
  // set prefix, file("peakOut/${prefix}.pval_signal.bw"), mark, fragLen, val("pvalueSignal") into peakCallWithInputResults

  script:
  //extSize = Math.round((fragLen as int)/2)
  """
  # narrow peaks and preliminary signal tracks
  macs2 callpeak -t ${bam} ${ control.exists() ? "-c ${control}" : "" } -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 --nomodel --extsize=${fragLen} \
                 --keep-dup all -B --SPMR
  # Broad and Gapped Peaks
  macs2 callpeak -t ${bam} ${ control.exists() ? "-c ${control}" : "" } -n ${prefix} --outdir peakOut \
                 -f BAM -g ${params.genomeSize} -p 1e-2 --broad --nomodel --extsize=${fragLen} \
                 --keep-dup all
  # rescale peaks on 10-1000 scale
  ${ params.rescale ? 
    ['narrow', 'broad', 'gapped'].collect { type ->
      def rescale_awk_str = 'BEGIN{FS=OFS="\\t";min=1e20;max=-1}'
      rescale_awk_str += 'NR==FNR&&NF!=0{min>$5?min=$5:min=min;max<$5?max=$5:max=max;next}'
      rescale_awk_str += 'NF!=0{n=$5;x=10;y=1000;$5=int(((n-min)*(y-x)/(max-min))+x);print}'
      "awk '${rescale_awk_str}' peakOut/${prefix}_peaks.${type}Peak peakOut/${prefix}_peaks.${type}Peak > peakOut/${prefix}_peaks.${type}Peak_rescaled && mv peakOut/${prefix}_peaks.${type}Peak{_rescaled,}"
    } : ""
  }
  # pileup signal tracks
  slopBed -i peakOut/${prefix}_treat_pileup.bdg -g ${chromSizes} -b 0 \
  | bedClip stdin ${chromSizes} peakOut/${prefix}.pileup_signal.bedgraph
  bedGraphToBigWig peakOut/${prefix}.pileup_signal.bedgraph ${chromSizes} peakOut/${prefix}.pileup_signal.bw
  ${ control.exists() ? 
    "# Fold enrichment signal tracks
    macs2 bdgcmp -t peakOut/${prefix}_treat_pileup.bdg \
  	             -c peakOut/${prefix}_control_lambda.bdg --outdir peakOut \
                 -o ${prefix}_FE.bdg -m FE
    slopBed -i peakOut/${prefix}_FE.bdg -g ${chromSizes} -b 0 \
     | bedClip stdin ${chromSizes} peakOut/${prefix}.fc.signal.bedgraph
    bedGraphToBigWig peakOut/${prefix}.fc.signal.bedgraph ${chromSizes} peakOut/${prefix}.fc_signal.bw
    # -log10(p-value) signal tracks
    bamReads=\$(samtools view -c ${bam}) && controlReads=\$(samtools view -c ${control}) \
     && sval=\$(bc -l <<< \$((bamReads<controlReads?bamReads:controlReads))/1000000)
    macs2 bdgcmp -t peakOut/${prefix}_treat_pileup.bdg \
  	             -c peakOut/${prefix}_control_lambda.bdg --outdir peakOut \
                 -o ${prefix}_ppois.bdg -m ppois -S \$sval
    slopBed -i peakOut/${prefix}_ppois.bdg -g ${chromSizes} -b 0 \
     | bedClip stdin ${chromSizes} peakOut/${prefix}.pval.signal.bedgraph
    bedGraphToBigWig peakOut/${prefix}.pval.signal.bedgraph ${chromSizes} peakOut/${prefix}.pval_signal.bw"
    : ""
  }
  rm -rf peakOut/${prefix}*.bdg peakOut/${prefix}*.bedgraph peakOut/${prefix}*.xls peakOut/${prefix}*.bed
  """
}

// process peakCallNoInput {
//   input:
//   file chromSizes from chromSizesForNoInput.val
//   set prefix, file(bam), mark, fragLen, view from bamsNoInput

//   output:
//   set prefix, file("peakOut/${prefix}_peaks.narrowPeak"), mark, fragLen, val("narrowPeak") into peakCallNoInputResults
//   set prefix, file("peakOut/${prefix}_peaks.broadPeak"), mark, fragLen, val("broadPeak") into peakCallNoInputResults
//   set prefix, file("peakOut/${prefix}_peaks.gappedPeak"), mark, fragLen, val("gappedPeak") into peakCallNoInputResults
//   set prefix, file("peakOut/${prefix}.pileup_signal.bw"), mark, fragLen, val("pileupSignal") into peakCallNoInputResults

//   script:
//   //extSize = Math.round((fragLen as int)/2)
//   command = ""
//   // narrow peaks and preliminary signal tracks
//   command += "macs2 callpeak -t ${bam} -n ${prefix} --outdir peakOut"
//   command += " -f BAM -g ${params.genomeSize} -p 1e-2 --nomodel --extsize=${fragLen}"
//   command += " --keep-dup all -B --SPMR\n"
//   // Broad and Gapped Peaks
//   command += "macs2 callpeak -t ${bam} -n ${prefix} --outdir peakOut"
//   command += " -f BAM -g ${params.genomeSize} -p 1e-2 --broad --nomodel --extsize=${fragLen}"
//   command += " --keep-dup all\n"
//   // rescale peaks on 10-1000 scale
//   if ( params.rescale ) {
//     ['narrow', 'broad', 'gapped'].collect { type ->
//       rescale_awk_str = 'BEGIN{FS=OFS="\\t";min=1e20;max=-1}'
//       rescale_awk_str += 'NR==FNR&&NF!=0{min>$5?min=$5:min=min;max<$5?max=$5:max=max;next}'
//       rescale_awk_str += 'NF!=0{n=$5;x=10;y=1000;$5=int(((n-min)*(y-x)/(max-min))+x);print}'
//       command += "awk '${rescale_awk_str}' peakOut/${prefix}_peaks.${type}Peak peakOut/${prefix}_peaks.${type}Peak"
//       command += " > peakOut/${prefix}_peaks.${type}Peak_rescaled && mv peakOut/${prefix}_peaks.${type}Peak{_rescaled,}\n"
//     }
//   }
//   // pileup signal tracks
//   command += "slopBed -i peakOut/${prefix}_treat_pileup.bdg -g ${chromSizes} -b 0"
//   command += " | bedClip stdin ${chromSizes} peakOut/${prefix}.pileup_signal.bedgraph\n"
//   command += "bedGraphToBigWig peakOut/${prefix}.pileup_signal.bedgraph ${chromSizes} peakOut/${prefix}.pileup_signal.bw\n"
//   command += "rm -rf peakOut/${prefix}*.bdg peakOut/${prefix}*.bedgraph peakOut/${prefix}*.xls peakOut/${prefix}*.bed"
// }

process NRF {
  input:
  set prefix, file(bam), controlId, mark, view from bams4NRF

  output:
  set prefix, stdout into NRFBams

  script:
  cpus = task.cpus
  """
  samtools view -@${cpus} ${bam} | NRF.awk 
  """
}

(peakCallResults, peakCallResults4FRiP) = peakCallWithInputResults
.flatMap { prefix, files, mark, fragLen ->
   def peaks = []
   files.each { f ->
        def view = (f.name =~ /([^\.]+Peak|[^.]+signal)/)[0][1].tokenize('_').indexed().collect { i, it->
          i == 0 ? it : it.toLowerCase().capitalize()
        }.join("")
        peaks << [ prefix, files, mark, fragLen, view ]
    }
    peaks
}.into(2)

input4FRiP = bams4FRiP.map { prefix, bam, controlId, mark, fragLen, view ->
    [prefix, bam, mark, fragLen, view]
}.mix(peakCallResults4FRiP.filter { it ->
  it[4] == 'narrowPeak'
}).groupTuple(by: [0,2,3])
.map { prefix, files, controlId, mark, views ->
  [prefix, files[0], files[1]]
}

process FRiP {

  input:
  set prefix, file(bam), file(peak) from input4FRiP

  output:
  set prefix, stdout into FRiPBams

  script:
  cpus = task.cpus
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
results.map { prefix, bam, control, mark, fragLen, view ->
  [ prefix, bam, mark, fragLen, view ]
}
.mix(peakCallResults))
.map { qc, result ->
    result + qc[1..-1]
}
.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) { prefix, path, mark, fragLen, view, nrf, frip ->
    [ prefix, path, mark, fragLen, view, nrf, frip ].join("\t")
}
.subscribe {
    log.info ""
    log.info "-----------------------"
    log.info "Pipeline run completed."
    log.info "-----------------------"
}
