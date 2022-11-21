/*
General workflow, Take contigs as input and output bgcprediction fasta/gbk
*/

//import modules
include {filterbysize}			from  "../modules/filter"
include {filter_pfamgbk}		from 	"../modules/filter"
include {splitbysize}				from	"../modules/filter"
include {deepbgc_prepare} 	from	"../modules/bgc"
include	{deepbgc_detect}		from  "../modules/bgc"
include {splitgbk}					from	"../modules/bgcgbk"
include {runsanntis}				from  "../modules/bgc"
include {runantismash}			from  "../modules/bgc"
include {rungecco}					from  "../modules/bgc"
include {prodigal}					from	"../modules/genepredict"
include {interproscan}			from	"../modules/genepredict"
include {sanntisgbk}				from	"../modules/genepredict"

//run metagenomic assembly pipeline using megahit
process modify_contigid_splitfas {
input:
  tuple val(x), path(splitcontig)
output:
  tuple val(split_name),  path(splitcontig), emit: splitcontigs
script:
  split_name=splitcontig.getName() 
"""
echo ${split_name}
"""
}



workflow BGCPRED {
	folder=params.contig_folder
	fasta=Channel.fromPath(["${folder}/*contigs.fa", "${folder}/*.fasta", "${folder}/*.fna"])
	filterbysize(fasta)
	longcontigs=filterbysize.out.contigs
	sizecontigs=longcontigs.splitFasta(size: "80.MB" ,file:true)
	modify_contigid_splitfas(sizecontigs)
	bysizecontigs=modify_contigid_splitfas.out.splitcontigs
	// Split large contig file into ~100 MB files
	// Bedfiles chomosome name, start, end, feature name

	if(params.antismash){
		runantismash(bysizecontigs)
		gbk_as=runantismash.out.gbk
		//as2bed(gbk_as)
		//bed_as=runantismash.out.bed
	}

	if(params.sanntis){
		prodigal(bysizecontigs)
		genesfaa=prodigal.out.genesfaa
		intogbk=bysizecontigs.join(genesfaa)
		sanntisgbk(intogbk)
		gbk=sanntisgbk.out.gbk
		interproscan(genesfaa)
		gff3=interproscan.out.iptsv
		sanntisinput=gff3.join(gbk)
		sanntisinput.view()
		runsanntis(sanntisinput)
		gff_sn=runsanntis.out.gff
		
		//runsanntis(bysizecontigs)
		//bed_sn=runsanntis.out.bed
		//fasta_sn=runsanntis.out.fasta //convert
		//gff2gbk(gff_sn, fasta_sn)
		//gbk_sn=gcf2gbk.out.gbk //convert into gbk and reformat for BigSCAPE
		//bedops workflow, receive bedfiles with gbk names, and genbank files
		// incremental addition of the next bgc predictors
		//bedops(bed_as, bed_sn)
	}

	if(params.gecco){
		//skash/gecco-0.6.3:latest
		rungecco(bysizecontigs) //remove predictions on edges
		gbk_g1=rungecco.out.gbk
		tsv_gc=rungecco.out.tsv
		//bgc_gc=
		//reformatgbk(gbk_g1)
		//gbk_gc=reformatgbk.out.gbk	//reformat for BiGSCAPE
	}

	if(params.deepbgc){
	  	deepbgc_prepare(bysizecontigs)
		contig_gbk=deepbgc_prepare.out.gbk
		filter_pfamgbk(contig_gbk)
		gbk=filter_pfamgbk.out.pf_gbk
 		deepbgc_detect(gbk)
		pred_bgcgbk=deepbgc_detect.out.bgc_gbk
		//reformat into Antishmash format
		//reformatgbk(pred_bgcgbk)
		//gbk_dp=regormatgbk.out.gbk
	}

	//transform into fasta, mash95, mcl
	// prioritize AntismashOutput, longest sequence_length
	// score different pipelines
	// Build database
	// Compare AS-sanntis -> AS(1), get intersection-predition tool table
	// Compare AS(1)-gecco -> AS(2), get intersection
	// Compare AS(2)-deepbgc -> AS(3), get intersection
	//The comparison steps: Overlap between BGC prediction
	// Options for finding the overlap:
	// compare final predicted sequecnes and keep:
	// 1 - all AS predictions
	// 2 - only intersections with other predictiors >1
	// 3 - only very high scoring predicions single predictor

	// convert gbk to fasta for each gbk in contig file,
	// work with list of individual files gbk or fasta.
	// mash distace each AS file with 2darybgc predictos
	// mask location of BGC prediction in contig file with
	// use bedtools create/identify features in contig file and
	// ask for intersectoin, return intersection features,

	////Find overlap using bed files, use bedops and bedfiles
	// Reference Genome are the contigs
	// The features are the GBK start and end, use this info to include
	// Create Bedfile from Predicted GBKs
	// info prediction on contig edge
	// Increase cutoff value for all predictors,
	// Create genome bed file for contigs file
	// P
	//antismash(fnafilt)
}
