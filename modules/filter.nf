
process filterbysize{
  //set directives to control
  //module "bioinfo-tools: SeqKit"
  scratch true
  container  "quay.io/biocontainers/seqkit:0.15.0--0"
  // memory '6GB'
  cpus '1'
  time '30m'

  input:
    path(fasta)
  output:
    tuple val(samp_name), path ("${samp_name}.fna") , emit: contigs

  script:
  samp_name=fasta.getSimpleName()

  """
  seqkit seq \\
  --min-len 5000  \\
  ${fasta} > ${samp_name}.filtzise.fna

  seqkit replace \\
  --pattern '(^)' \\
  --replacement '${samp_name}_' \\
  ${samp_name}.filtzise.fna > ${samp_name}.fna

  """
}

//the module filters out gbk files with no predicted pfams into it
process filter_pfamgbk{
  //set directives
  container "biopython/biopython:latest"
  //module "python3"
  //conda "conda-forge::python=3.6.7 conda-forge::biopython=1.74"
  errorStrategy {task.exitStatus in 1 ? 'ignore': 'retry'}
  maxRetries = 3
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

//The module takes a contig file, checks the file size and divides
//it by 70MB, returns the estimated number of files.
process splitbysize{
  //set directives
  container 'quay.io/biocontainers/pyfasta:0.5.2--py_1'
  scratch true
  cpus '1'
  time '15m'

  input:
  tuple val(x), path(contigs)
  output:
  tuple val(x),path("*.fasta"), emit: contigs

  script:
  """
  filename=${contigs}
  filesize=\$(stat -c%s "\$filename")
  float=calc \$filesize/10000000
  int=\${float%.*}
  if ((int > 1));
    then
    pyfasta split -n \$int ${contigs}
    else
    cat ${contigs} > ${x}.1.fasta
  fi
  """
}
