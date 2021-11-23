workflow parseIndexFile {
  take:
    index
    globalFragLen
  main:
    log.info("Parsing index file")
    fastqs = Channel
      .from(index.readLines())
      .map { line ->
          def list = line.tokenize()
          def mergeId = list[0]
          def id = list[1]
          def path = resolveFile(list[2], index)
          def controlId = list[3]
          def mark = list[4]
          def fragLen = params.fragmentLength
          if ( params.shift || !fragLen ) {
              fragLen = list[5] as Integer
          }
          def quality = fastq(path).qualityScore()
          [ mergeId, id, path, controlId, mark, fragLen, quality ]
      }
      printFragInfo(fastqs, globalFragLen)
    emit:
      fastqs
}

def printFragInfo(fastqs, globalFragLen) {
  colWidth = fastqs.map {
      it[0].size()
    }.max()
    fragInfo = fastqs
      .combine(colWidth)
      .map {
      [ it[0], params.shift, it[5], it[-1] ]
    }.view { mergeId, global, fragLen, colWidth ->
      def message = "[${mergeId.padRight(colWidth)}] > "
      if ( global ) {
          message += "global fragment length: ${globalFragLen} - shift size by: "
      } else {
        message += "fragment length: "
      }
      if ( fragLen ) {
          message += "${fragLen}"
      } else {
          message += "estimated"
      }
      message
    }
}

/*
 * Given a string path resolve it against the index file location.
 * Params:
 * - str: a string value represting the file path to be resolved
 * - index: path location against which relative paths need to be resolved
 */
def resolveFile( str, index ) {
  if( str.startsWith('/') || str =~ /^[\w\d]*:\// ) {
    return file(str)
  }
  else if( index instanceof Path ) {
    return index.parent.resolve(str)
  }
  else {
    return file(str)
  }
}


def printUsage(defaults) {
    //print usage
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
    log.info '    --genome-index GENOME_INDEX_FILE    Reference genome index file.'
    log.info '    --genome-size GENOME_SIZE           Reference genome size for MACS2 callpeaks. Must be one of'
    log.info "                                        MACS2 precomputed sizes: hs, mm, dm, ce. (Default: '${defaults.genomeSize}')"
    log.info "    --db-file FILE                      Output file where to store results information (Default: '${defaults.dbFile}')"
    log.info "    --replicate-pattern PATTERN         Glob pattern used to match replicates (Default: '${defaults.replicatePattern}')."
    log.info "    --mismatches MISMATCHES             Sets the maximum number/percentage of mismatches allowed for a read (Default: '${defaults.mismatches}').  "
    log.info "    --multimaps MULTIMAPS               Sets the maximum number of mappings allowed for a read (Default: '${defaults.multimaps}')."
    log.info "    --min-matched-bases BASES           Sets the minimum number/percentage of bases that have to match with the reference (Default: '${defaults.minMatchedBases}')."
    log.info "    --quality-threshold THRESHOLD       Sets the sequence quality threshold for a base to be considered as low-quality (Default: '${defaults.qualityThreshold}')."
    log.info "    --fragment-length LENGTH            Sets the fragment length globally for all samples (Default: '${defaults.fragmentLength}')."
    log.info "    --zerone-min-confidence CONFIDENCE  Make Zerone print targets with confidence higher than CONFIDENCE (Default: ${defaults.zeroneMinConfidence})."
    log.info "    --remove-duplicates                 Remove duplicate alignments instead of just flagging them (Default: '${defaults.removeDuplicates}')."
    log.info '    --rescale                           Rescale peak scores to conform to the format supported by the'
    log.info "                                        UCSC genome browser (score must be <1000) (Default: '${defaults.rescale}')."
    log.info "    --shift                             Move fragments ends and apply global extsize in peak calling (Default: '${defaults.shift}')."
    log.info "                                        If '--shift' is set and '--fragment-length' is not sepcified the global fragmenth length"
    log.info "                                        is forced to '200'."
    log.info ''
}

def printHeader() {
    ////// Print parameters ///////
    log.info ''
    log.info 'C H I P - N F ~ ChIP-seq Pipeline'
    log.info '---------------------------------'
    log.info ''
    log.info "General parameters"
    log.info '------------------'
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
    log.info "Max Mismatches         : ${params.mismatches}"
    log.info "Max Multimaps          : ${params.multimaps}"
    log.info "Minimum Matched Bases  : ${params.minMatchedBases}"
    log.info "Low Quality Threshold  : ${params.qualityThreshold}"
    log.info ''
    log.info "Zerone parameters"
    log.info '-----------------'
    log.info "Confidence Threshold   : ${params.zeroneMinConfidence}"
    log.info ''
}

// Collect results data and write it into DB
workflow collectResults {
    take:
      metrics
      bams
      peaks
      signal
      zerone
    main:
      bamResults = bams.map {
        it[0..1] + it[3..-2]
      }
      zeroneResults = zerone.map {
        it[0..2] + ['-', it[3], '-', '-']
      }
      inputResults = bamResults.mix(signal)
        .filter { it[2] == 'input'}
        .map {
          it + ['-', '-']
        }
      results = metrics.cross(
        bamResults.mix(peaks).mix(signal)
      ).map { qc, result ->
        result + qc[1..-1]
      }
      .mix(inputResults)
      .mix(zeroneResults)

      //Init DB
      pdb = file(params.dbFile)
      pdb.write('')

      // Collect results
      results.collectFile(name: pdb.name, storeDir: pdb.parent, newLine: true) {
        // prefix, path, mark, fragLen, view, nrf, frip
        it.join("\t")
      }
}
