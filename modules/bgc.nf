process deepbgc_prepare{
  scratch=true
  cpus 2
  container 'quay.io/biocontainers/deepbgc:0.1.27--pyhdfd78af_0'
  publishDir "deepbgc/prepare"
  time   = { 30.h  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(fasta)
  output:
    tuple val(x), path("${x}.csv"), emit: pfamcsv, optional: true
<<<<<<< HEAD
    tuple val(x), path("${x}_full.gbk" ), emit: gbk, optional: true
=======
    tuple val(x), path("${x}_full.gbk"), emit: gbk, optional: true
>>>>>>> 8b1d6080c0b97f8e668c9f10affa1b3b41ece644
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
  errorStrategy {task.exitStatus in 1 ? 'ignore': 'terminate'}
  //validExitStatus 0,1

  input:
    tuple val(x), path(pfamgbk)
  output:
<<<<<<< HEAD
    tuple val(x), path("${x}/*.bgc.gbk"), emit: bgc_gbk, optional:true
    tuple val(x), path("${x}/*.bgc.tsv"), emit: bgc_tsv optional:true
=======
    tuple val(x), path("${x}/*.bgc.gbk"), optional:true, emit: bgcgbk
    path("${x}/*.bgc.tsv"), optional:true, emit: tsv
    path("${x}/*.json"), optional:true,   emit: json
>>>>>>> 8b1d6080c0b97f8e668c9f10affa1b3b41ece644
  script:

    """
    deepbgc pipeline \\
    	--output ${x} \\
<<<<<<< HEAD
      --label deepbgc_80_score \\
      --score 0.8 \\
      --min-proteins 2 \\
      ${pfamgbk}
=======
      	--label deepbgc_80_score \\
      	--score 0.8 \\
      	${pfamgbk}
>>>>>>> 8b1d6080c0b97f8e668c9f10affa1b3b41ece644
    """
}
