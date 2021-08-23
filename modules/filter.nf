
process filterbysize{
  //set directives to control
  scratch true
  memory '6GB'
  cpus '1'
  time '30m'
  module 'bioinfo-tools:SeqKit'

  input:
    path(fasta)
  output:
    tuple val(samp_name), path ("${samp_name}.filt.contig.fna") , emit: contigs

  script:
  samp_name=fasta.getSimpleName()

  """
  seqkit seq --min-len 5000  ${fasta} > ${samp_name}.filt.contig.fna

  """
}
