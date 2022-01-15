process splitgbk{
  scratch=true
  cpus 1
  time '10m'
  publishDir "deepbgc_split"

  input:
    tuple val(x), path(gbk)
  output:
    tuple val(x), path("*.gbk" ), emit: bgcgbk
  script:
    """
    splitgbk $gbk
    rm $gbk
    """
}
