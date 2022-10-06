process fastaIndex {
  input:
  file genome

  output:
  file "${genomePrefix}.fai"

  script:
  def cat = genome.name.endsWith('.gz') ? 'zcat' : 'cat'
  genomePrefix = genome.name.endsWith('.gz') ? genome.baseName : genome.name
  """
  mkfifo genome_fifo.fa
  ${cat} ${genome} > genome_fifo.fa &
  samtools faidx genome_fifo.fa -o ${genomePrefix}.fai
  """
}
