/*
General workflow, Take contigs as input and output bgcprediction fasta/gbk
*/

//import modules
<<<<<<< HEAD
include {filterbysize}		from  "../modules/filter"
include {deepbgc_prepare}	from  "../modules/bgc"
include {deepbgc_detect}	from  "../modules/bgc"
=======
include {filterbysize}			from  "../modules/filter"
include {deepbgc_prepare} 	from	"../modules/bgc"
include	{deepbgc_detect}		from  "../modules/bgc"
>>>>>>> 265d9cb040bd1c7fe707ca424bb0e256acf4c263
//include {antismash}		from  "../modules/bgc"

//run metagenomic assembly pipeline using megahit

workflow BGCPRED {
	folder=params.contig_folder
	fasta=Channel.fromPath(["${folder}/*contigs.fa", "${folder}/*.fasta", "${folder}/*.fna"])
	filterbysize(fasta)
	longcontigs=filterbysize.out.contigs
  deepbgc_prepare(longcontigs)
	contig_pfamgbk=deepbgc.out.gbk
	deepbgc_detect(contig_pfamgbk)

	//antismash(fnafilt)
        }
