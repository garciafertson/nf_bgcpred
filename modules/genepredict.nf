process prodigal {
  cpus '1'
  time '4h'
  container 'nanozoo/prodigal:2.6.3--2769024'

  input:
    tuple val(samp_name), path(contig)
  output:
    tuple val(samp_name), path("*.prodigal.faa"), emit: genesfaa

  script:
    //samp_name=contig.getSimpleName()
    """
    prodigal -i ${contig} \\
    -a ${samp_name}.prodigal.faa \\
    -p meta
    """
}

process interproscan{
  cpus '2'
  time '10h'
  container 'interpro/interproscan:5.51-85.0'

  input:
    tuple val(x), path(genesfaa)
  output:
    tuple val(x), path("*.ip.gff3"), emit: ipgff3
  script:
  """
  interproscan.sh -i ${genesfaa} \\
  -f GFF3 \\
  -appl Pfam,TIGRFAM,PRINTS,ProSitePatterns,Gene3D \\
  -cpu 2 \\
  -o ${x.id}.ip.gff3
  """
}

process sanntisgbk{
  cpus '2'
  time '10h'
  container 'sysbiojfgg/sanntis:0.1'

  input:
  tuple val(x), path(fna), path(genesfaa)
  output:
  tuple val(x), fpath("*.faa.gb"), emit: genebank
  script:
  """
  sanntis_build_gb -n ${fna}
    -a ${genesfaa}
    -o ${x.id}.faa.gb
  """
}
