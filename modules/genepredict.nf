process prodigal {
  cpus '1'
  time '4h'
  container 'nanozoo/prodigal:2.6.3--2769024'

  input:
    tuple val(samp_name), path(contig)
  output:
    tuple val(samp_name), path("*.prodigal.faa"), emit: genesfaa

  script:
    """
    prodigal -i ${contig} \\
    -a ${samp_name}.faa \\
    -p metai
    sed 's/*//' ${samp_name}.faa > ${samp_name}.prodigal.faa
    """
}

process interproscan{
  cpus '3'
  time '15h'
  container 'interpro/interproscan:5.51-85.0'

  input:
    tuple val(x), path(genesfaa)
  output:
    //tuple val(x), path("*.ip.gff3"), emit: ipgff3
    tuple val(x), path("*ip.tsv"), emit: iptsv
  script:
  """
  ls /opt/interproscan/data
  interproscan.sh -i ${genesfaa} \\
  -f TSV \\
  -appl Pfam,TIGRFAM,PRINTS,ProSitePatterns,Gene3D \\
  -cpu 6 \\
  -o ${x}.ip.tsv
  """
}

process sanntisgbk{
  cpus '1'
  time '1h'
  container 'sysbiojfgg/sanntis:0.1'

  input:
  tuple val(x), path(fna), path(genesfaa)
  output:
  tuple val(x), path("*.faa.gb"), emit: gbk
  script:
  """
  sanntis_build_gb -n ${fna} \\
    -a ${genesfaa} \\
    -o ${x}.faa.gb
  """
}
