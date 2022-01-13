/*
General workflow, Take contigs as input and output bgcprediction fasta/gbk
*/

//import modules
include {filterbysize}		from  "../modules/filter"
include {deepbgc_prepare}	from  "../modules/bgc"
include {deepbgc_detect}	from  "../modules/bgc"
//include {antismash}		from  "../modules/bgc"

//run metagenomic assembly pipeline using megahit

workflow BGCPRED {
	folder=params.contig_folder
	fasta=Channel.fromPath(["${folder}/*contigs.fa", "${folder}/*.fasta", "${folder}/*.fna"])
	filterbysize(fasta)
	longcontigs=filterbysize.out.contigs
  	deepbgc_prepare(longcontigs)
	contig_pfamgbk=deepbgc_prepare.out.gbk
	deepbgc_detect(contig_pfamgbk)
	//antismash(fnafilt)
        }
