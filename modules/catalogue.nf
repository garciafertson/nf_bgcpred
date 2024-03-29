process bedops_dp{
  //scratch=true
  cpus 1
  container 'biocontainers/bedops:v2.4.35dfsg-1-deb_cv1'
  time   = { 20.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(bed_ref), path(bed_add)
  output:
    tuple val(x), path("${x}_jdp.bed"), emit: bed
  script:
    """
    sort ${bed_add} > add
    sort ${bed_ref} > ref
    bedops --not-element-of 20% add ref > ${x}new.bed
    cat ${bed_ref} ${x}new.bed > ${x}_jdp.bed
    """
}

process bedops_gc{
  //scratch=true
  cpus 1
  container 'biocontainers/bedops:v2.4.35dfsg-1-deb_cv1'
  time   = { 20.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(bed_ref), path(bed_add)
  output:
    tuple val(x), path("${x}_jgc.bed"), emit: bed
  script:
    """
    sort ${bed_add} > add
    sort ${bed_ref} > ref
    bedops --not-element-of 20% add ref > ${x}new.bed
    cat ${bed_ref} ${x}new.bed > ${x}_jgc.bed
    """
}

process bedops_sn{
  //scratch=true
  cpus 1
  container 'biocontainers/bedops:v2.4.35dfsg-1-deb_cv1'
  time   = { 20.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(bed_ref), path(bed_add)
  output:
    tuple val(x), path("${x}_jsn.bed"), emit: bed
  script:
    """
    sort ${bed_add} > add
    sort ${bed_ref} > ref
    bedops --not-element-of 20% add ref > ${x}new.bed
    cat ${bed_ref} ${x}new.bed > ${x}_jsn.bed
    """
}

process getbgc_fna{
  //scratch=true
  cpus 1
  container 'staphb/bedtools:2.30.0'
  publishDir "out/allbgc", mode:"copy", overwrite:true
  time   = { 40.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2

  input:
    tuple val(x), path(bgcbed), path(contigs)
  output:
    path("*_bgc.fna"), emit: fna
  script:
    """
    bedtools getfasta -fi ${contigs} -bed ${bgcbed} -nameOnly > ${x}_bgc.fna
    """
}

process getbgc_bed{
  //scratch=true
  cpus 1
  time   = { 40.m  * task.attempt }
  errorStrategy = 'retry'
  maxRetries = 2
  publishDir "out/allbgc", mode:"copy", overwrite:true

  input:
    tuple val(x), path(bgcbed)
  output:
    path("*_bgc.bed"), emit: bed
  script:
    """
    cat ${bgcbed} > ${x}_bgc.bed
    """
}


process mashtriangle{
  //directives
  container "staphb/mash:2.3"
  cpus 10
  time 24.h

  input:
    path(mashfile)
    output:
    path("*.edgelist"), emit: edgelist
    path("*_0.05.list"), emit: list05
  script:
      """
      mash triangle -p 10 -E -d 0.3 -i $mashfile > ${mashfile}.edgelist
      awk '{if (\$3 < $params.mashdistance) print \$1,\$2,\$3}' \\
      ${mashfile}.edgelist > ${mashfile}_0.05.list
      """
}

process mcl_clust{
  //directives
  //module "mcl"
  container "sysbiojfgg/mcl:v0.1"
  publishDir "out/bgc_catalogue", mode:'copy', overwrite:true
  cpus 4
  time 5.h

  input:
    path(distances)
  output:
    path("*.clusters"), optional: true, emit: clusters
  script:
    """
    if [ -s $distances ]; then
    mcl $distances --abc -te 4 -I $params.inflation -o ${distances}.clusters;
    else
    touch ${distances}.clusters
    fi
    """
}

process fna_get_representatives{
  //directives
  publishDir "out/bgc_catalogue", mode: 'copy', overwrite:true
  container "biopython/biopython:latest"
  cpus 1
  time 4.h

  input:
    path(bed)
    path(clusters)
    path(fna)
  output:
    path("*representative.bed"), emit: representative_bed
    path("*representative.fna"), emit: representative_fna
  script:
    """
    get_representatives.py  --bedfile ${bed} \\
                --clusters ${clusters} \\
                --fna ${fna} \\
    """
}

process build_index{
  //directives for uppmax rackham
  publishDir "out/bgc_catalogue",
  mode: "copy", overwrite:true
  //conda "bioconda::bowtie2"
  container "biocontainers/bowtie2:v2.4.1_cv1"
  cpus params.mashcores
  time 6.h

  input:
    path(seqs)

  output:
    path("bgc_catalogue_bowtieindex*"), emit: index

  script:
    """
    bowtie2-build -f --threads $task.cpus $seqs bgc_catalogue_bowtieindex
    """
}

process build_genomebed{
	//directives
	publishDir "out/bgc_catalogue",
		mode: "copy", overwrite:true
	cpus 1
	time 1.h
	//module "python3"
	container "biopython/biopython:latest"

	input:
	path(bed)
	output:
	path("bgc_catalogue.genome.bed"), emit: gnmbed
	script:
		"""
		genomebed.py $bed
		"""
}

process prodigal_bgc{
  cpus 1
  time 4.h
  //module "python3"
  container 'nanozoo/prodigal:2.6.3--2769024'

  input:
  path(fna)
  output:
  path("representative_bgc.gff"), emit: gff
  script:
    """
    prodigal \\
    -i ${fna} \\
    -c -f gff \\
    -o representative_bgc.gff \\
    -p meta
    """
}

process build_bedbgc{
  publishDir "out/bgc_catalogue",
    mode: "copy", overwrite:true
  cpus 1
  time 4.h
  //module "python3"
  container "biopython/biopython:latest"

  input:
  path(gff)
  output:
  path("bgc_catalogue.bed"), emit: gff
  script:
  """
  gff2bed.py \\
  --gff ${gff} \\
  --name bgc_catalogue
  """
}
