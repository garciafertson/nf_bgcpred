process deepbgc_prepare{
  scratch=true
  cpus 2
  container 'quay.io/biocontainers/deepbgc:0.1.27--pyhdfd78af_0'
  publishDir "deepbgc/prepare"
  //time   = { 35.h  * task.attempt }
  time=36.h
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(fasta)
  output:
    tuple val(x), path("${x}.csv"), emit: pfamcsv, optional: true
    tuple val(x), path("${x}_full.gbk"), emit: gbk, optional: true
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
  publishDir "deepbgc",
    mode: "copy",
    overwrite: true
  errorStrategy {task.exitStatus in 1 ? 'ignore': 'retry'}
  maxRetries = 3
  //validExitStatus 0,1

  input:
    tuple val(x), path(pfamgbk)
  output:
    tuple val(x), path("${x}/*.bgc.gbk"), optional:true, emit: bgc_gbk
    path("${x}/*.bgc.tsv"), optional:true, emit: tsv
    path("${x}/*.json"), optional:true,   emit: json
  script:

    """
    deepbgc pipeline \\
    	--output ${x} \\
      --label deepbgc_80_score \\
      --score 0.8 \\
      --min-proteins 2 \\
      ${pfamgbk}
    """
}
