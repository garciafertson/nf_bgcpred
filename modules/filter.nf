
process filterbysize{
  //set directives to control
  module "bioinfo-tools: SeqKit"
  scratch true
  // memory '6GB'
  cpus '1'
  time '30m'

  input:
    path(fasta)
  output:
    tuple val(samp_name), path ("${samp_name}.filt.contig.fna") , emit: contigs

  script:
  samp_name=fasta.getSimpleName()

  """
  seqkit seq \\
  --min-len 5000  \\
  ${fasta} > ${samp_name}.filtzise.fna

  seqkit replace \\
  --pattern '(^)' \\
  --replacement '${samp_name}_' \\
  ${samp_name}.filtzise.fna > ${samp_name}.filt.contig.fna

  """
}

process filter_pfamgbk{
  //set directives
  scratch true
  memory'1G'
  cpus '1'
  time '10m'

  input:
  tuple val(x), path(gbk)
  output:
  tuple val(x), path("${x}.full_pf.gbk"), emit: pf_gbk, optional: true

  script:
  """
  pfamfilt_gbk.py ${gbk} ${x}.full_pf.gbk
  """
}
