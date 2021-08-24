process deepbgc_prepare{
  scratch=true
  cpus 2
  time '12h'
  container 'quay.io/biocontainers/deepbgc:0.1.27--pyhdfd78af_0'
  publishDir "bgcprepare"

  input:
    tuple val(x), path(fasta)
  output:
    tuple val(x), path("${x}.csv"), emit: pfamcsv
    tuple val(x), path("${x}_full.gbk"), emit: gbk
  script:
    """
    deepbgc prepare \\
    --prodigal-meta-mode \\
    --output-gbk ${x}_full.gbk \\
    --output-tsv ${x}.csv \\
    ${fasta}
    """
}

process deepbgc_detect{
  scratch true
  cpus 1
  time '2h'
  container 'quay.io/biocontainers/deepbgc:0.1.27--pyhdfd78af_0'
  publishDir "deepbgc"

  input:
    tuple val(x), path(pfamgbk)
  output:
    tuple val(x), path("${x}/*.gbk"), emit: bgc_gbk
  script:

    """
    deepbgc pipeline \\
    	--minimal-output \\
    	--output ${x} \\
      --label deepbgc_80_score \\
      --score 0.8 \\
      ${pfamgbk}
    """
}
