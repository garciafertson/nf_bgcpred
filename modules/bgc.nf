process deepbgc{
  scratch true
  cpus "8"
  time '6h'
  container 'quay.io/biocontainers/deepbgc:0.1.27--pyhdfd78af_0'
  publishDir "deepbgc"

  input:
    tuple val(x), path(fasta)
  output:
    tuple val(x), path("mySequence"), emit: reads
  script:

    """
    deepbgc pipeline \\
    	--minimal-output \\
    	--output ${x} \\
      --label deepbgc_90_score \\
      --score 0.9 \\
      ${fasta}

    """
}
