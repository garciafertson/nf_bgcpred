process parse_asresults{
  //scratch=true
  cpus 1
  container "biopython/biopython"
  time   = { 20.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(as_gbk)
  output:
  tuple val(x), path("*_as.bed"), emit: bed
  tuple val(x), path("*_as.fna"), emit: fna, optional: true
  script:
    """
    antismash_gbk_reformat.py \\
    --gbk $as_gbk \\
    --name ${x}
    """
}

process parse_snresults{
  //scratch=true
  cpus 1
  container "biopython/biopython"
  time   = { 20.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(gff_sn), path(contigs)
  output:
  tuple val(x), path("*_sn.bed"), emit: bed
  tuple val(x), path("*_sn.fna"), emit: fna, optional: true
  script:
    """
    sanntis_gbk_reformat.py \\
    --gff ${gff_sn} \\
    --contigs ${contigs} \\
    --name ${x}
    """
}

process parse_gcresults{
  //scratch=true
  cpus 1
  container "biopython/biopython"
  time   = { 20.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(gbk), path(tsv)
  output:
    tuple val(x), path("*_gc.bed"), emit: bed
    tuple val(x), path("*_gc.fna"), emit: fna, optional: true
  script:
    """
    gecco_gbk_reformat.py \\
    --gbk ${gbk} \\
    --tsv ${tsv} \\
    --name ${x} 
    """
}

process parse_dpresults{
  //scratch=true
  cpus 1
  container "biopython/biopython"
  time   = { 20.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(dp_tsv), path(dp_gbk)
  output:
    tuple val(x), path("*_dp.bed"), emit: bed
    tuple val(x), path("*_dp.fna"), emit: fna, optional: true
  script:
    """
    deepbgc_gbk_reformat.py \\
    --tsv $dp_tsv \\
    --gbk $dp_gbk \\
    --name $x\\
    """
}
