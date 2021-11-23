process indexing {
    input:
    file genome

    output:
    file "genome_index.gem"

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

process mapping {
  input:
  file index
  tuple val(mergeId), val(prefix), file(fastq), val(controlId), val(mark), val(fragLen), val(quality)

  output:
  tuple val(mergeId), val(prefix), file("${prefix}_primary.bam"), val(controlId), val(mark), val(fragLen), val('Alignments')

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
                  -o ${prefix}_primary.bam
  """
}

workflow align {
    take:
        genome
        fastqs
    main:
        index = getGenomeIndex(genome)
        mapping(index, fastqs)
    emit:
        mapping.out
}

def getGenomeIndex(genome) {
    if (params.genomeIndex)
        return file(params.genomeIndex)
    indexing(file(genome))
}
