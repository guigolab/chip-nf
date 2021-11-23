process fastaIndex {
  input:
  file genome

  output:
  file "${genome}.fai"

  script:
  """
  samtools faidx ${genome}
  """
}
