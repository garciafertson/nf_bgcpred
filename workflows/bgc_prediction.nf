/*
General workflow, Take contigs as input and output bgcprediction fasta/gbk
*/

//import modules
include {filterbysize}			from  "../modules/filter"
include {deepbgc_prepare} 	from	"../modules/bgc"
include	{deepbgc_detect}		from  "../modules/bgc"
include {splitgbk}					from	"../modules/bgcgbk"
//include {antismash}		from  "../modules/bgc"

//run metagenomic assembly pipeline using megahit

workflow BGCPRED {
	folder=params.contig_folder
	fasta=Channel.fromPath(["${folder}/*contigs.fa", "${folder}/*.fasta", "${folder}/*.fna"])
	filterbysize(fasta)
	longcontigs=filterbysize.out.contigs
  deepbgc_prepare(longcontigs)
	contig_gbk=deepbgc_prepare.out.gbk
	filter_pfamgbk(contig_gbk)
	gbk=filter_pfamgbk.out.pf_gbk
	deepbgc_detect(gbk)
	pred_bgcgbk=deepbgc_detect.out.bgc_gbk
	//splitgbk(pred_bgcgbk)

	//antismash(fnafilt)
        }
