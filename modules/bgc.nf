process deepbgc_prepare{
  scratch=true
  cpus 2
  container 'quay.io/biocontainers/deepbgc:0.1.30--pyhca03a8a_2'
  time   = { 20.h  * task.attempt }
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
  //scratch true
  cpus 1
  time '2h'
  container 'quay.io/biocontainers/deepbgc:0.1.30--pyhca03a8a_2'
  publishDir "out/deepbgc"
  errorStrategy {task.exitStatus in 1 ? 'ignore': 'terminate'}
  //validExitStatus 0,1

  input:
    tuple val(x), path(pfamgbk)
  output:
    tuple val(x), path("${x}_dp.gbk"), optional:true, emit: bgc_gbk
    path("${x}_dp.tsv"), optional:true, emit: tsv
  script:
    """
    deepbgc pipeline \\
    	--output ${x} \\
      --label deepbgc \\
      --score 0.925 \\
      --min-proteins 4 \\
      ${pfamgbk}
    mv ${x}/*.bgc.gbk ${x}_dp.gbk
    mv ${x}/*.bgc.tsv ${x}_dp.tsv
    """
}

process runantismash {
  //scratch true
  cpus 2
  time '10h'
  container 'antismash/standalone-nonfree:6.1.1'
  publishDir "out/antismash",
    mode:"copy",
    overwrite: true
  errorStrategy {task.exitStatus in 1 ? 'ignore': 'terminate'}
  //validExitStatus 0,1

  input:
  tuple val(x), path(contigs)
  output:
  tuple val(x), path("${x}/*.region???.gbk"), optional: true, emit: bgcgbk
  tuple val(x), path("${x}_as.gbk"), optional: true, emit:gbk
  script:
  """
  antismash --cb-general \\
  --cb-knownclusters \\
  --asf \\
  --genefinding-tool prodigal-m \\
  ${contigs}
  cat ${x}/*.region???.gbk > ${x}_as.gbk
  # rename
  """
}

process rungecco {
    //scratch true
    cpus 2
    time '10h'
    container 'skash/gecco-0.6.3:latest'
    errorStrategy {task.exitStatus in 1 ? 'ignore': 'terminate'}
    publishDir "out/gecco"
    //validExitStatus 0,1

    input:
    tuple val(x), path(contigs)
    output:
    tuple val(x), path("${x}.gc.gbk"), optional: true, emit: gbk
    tuple val(x), path("*clusters.tsv"), optional: true, emit: tsv
    //recover and analize clusters.tsv snd create bed in relation to conitgs(genome) file

    script:
    """
    gecco run --genome ${contigs} \\
    -o gecco_out \\
    --jobs $task.cpus \\
    --threshold 0.9
    mv gecco_out/*clusters.tsv ${x}_clusters.tsv
    cat gecco_out/*gbk > ${x}.gc.gbk
    """
}


process runsanntis {
    //scratch true
    cpus 1
    time '2h'
    container 'sysbiojfgg/sanntis:0.1'
    errorStrategy {task.exitStatus in 1 ? 'ignore': 'terminate'}
    //valifExitStatus 0,1
    //publishDir "out/sanntis"

    input:
    tuple val(x), path(tsv), path(gbk)
    output:
    tuple val(x), path("${x}sanntis/*ful.gff"), optional: true, emit: gff
    //recover and analize clusters.tsv snd create bed in relation to conitgs(genome) file

    """
    sanntis  \\
    --ip-file ${tsv} \\
    --antismash_output True \\
    ${gbk}
    """
    }
