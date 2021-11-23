// Enable DSL2
nextflow.enable.dsl=2

// Define default values for parameters
def defaults = [
  shift: false,
  mismatches: 2,
  multimaps: 10,
  rescale: false,
  genomeSize: 'hs',
  fragmentLength: 200,
  minMatchedBases: 0.8,
  qualityThreshold: 26,
  zeroneMinConfidence: 0,
  removeDuplicates: false,
  replicatePattern: '.[12]',
  dbFile: 'chipseq-pipeline.db',
  tmpDirMACS2: '.',
]

// Parameters for test run
params.index = "${baseDir}/data/index.tsv"
params.genome = "${baseDir}/data/genome.fa"

// Set default values for parameters
params.help = false
params.genomeIndex = ''
params.shift = defaults.shift
params.dbFile = defaults.dbFile
params.rescale = defaults.rescale
params.multimaps = defaults.multimaps
params.genomeSize = defaults.genomeSize
params.mismatches = defaults.mismatches
params.fragmentLength = 0
params.minMatchedBases = defaults.minMatchedBases
params.qualityThreshold = defaults.qualityThreshold
params.removeDuplicates = defaults.removeDuplicates
params.replicatePattern = defaults.replicatePattern
params.zeroneMinConfidence  = defaults.zeroneMinConfidence
params.tmpDirMACS2 = defaults.tmpDirMACS2

// Check input parameters
include { printUsage } from './modules/common'
if (params.help) {
  printUsage(defaults)
  exit 1
}
if (!params.genome) {
  exit 1, "Please specify a genome file"
}

if (!params.index) {
  exit 1, "Please specify the input table file"
}

// Includes
include { parseIndexFile; collectResults; printHeader } from './modules/common'
include { align } from './modules/gem'
include { mergeBams } from './modules/merge'
include { processBams } from './modules/bam'
include { computeMetrics } from './modules/metrics'
include { callPeaks as MACS2callPeaks} from './modules/macs'
include { callPeaks as ZERONEcallPeaks} from './modules/zerone'

workflow {
  printHeader()
  fastqs = parseIndexFile(
    file(params.index),
    params.fragmentLength ?: defaults.fragmentLength,
  )
  align(params.genome, fastqs)
  mergeBams(align.out)
  allBams = processBams(mergeBams.out)
  (peaks, signal) = MACS2callPeaks(
    file(params.genome),
    params.fragmentLength ?: defaults.fragmentLength,
    allBams
  )
  (zeroneMatrix, zeroneBed, zeroneMergedBed) = ZERONEcallPeaks(allBams)
  metrics = computeMetrics(
    allBams,
    peaks.filter { it[4] =~ /narrowPeak$/ }
  )
  collectResults(
    metrics,
    allBams,
    peaks,
    signal,
    zeroneMatrix.mix(zeroneBed).mix(zeroneMergedBed)
  )
}

workflow.onComplete {
    log.info ""
    log.info "-----------------------"
    log.info "Pipeline run completed."
    log.info "-----------------------"
}
