/*
General workflow, Take contigs as input and output bgcprediction fasta/gbk
*/

//import modules
include {filterbysize}			from  "../modules/filter"
include {filter_pfamgbk}		from 	"../modules/filter"
include {deepbgc_prepare} 	from	"../modules/bgc"
include	{deepbgc_detect}		from  "../modules/bgc"
include {splitgbk}					from	"../modules/bgcgbk"
include {runsanntis}				from  "../modules/bgc"
include {runantismash}			from  "../modules/bgc"
include {rungecco}					from  "../modules/bgc"

//run metagenomic assembly pipeline using megahit

workflow BGCPRED {
	folder=params.contig_folder
	fasta=Channel.fromPath(["${folder}/*contigs.fa", "${folder}/*.fasta", "${folder}/*.fna"])
	filterbysize(fasta)
	longcontigs=filterbysize.out.contigs

	// Derreplicate contig sequences MASH,
	derrep(longcontigs)
	derrepcontigs= derrep.out.contigs

	// Split large contig files into many files ~40MB
	splitbysize(derrepcontigs)
	bysizecontigs=splitbysize.out.contigs
	genome_bed(bysizecontigs)
	// Bedfiles chomosome name, start, end, feature name

	runantismash(bysizecontigs)
	gbk_as=runantismash.out.gbk
	//as2bed(gbk_as)
	//bed_as=runantismash.out.bed

	if(params.sanntis){
		runsanntis(bysizecontigs) //remove predictions on edges
		gff_sn=runsanntis.out.gff
		//bed_sn=runsanntis.out.bed
		//fasta_sn=runsanntis.out.fasta //convert
		//gff2gbk(gff_sn, fasta_sn)
		//gbk_sn=gcf2gbk.out.gbk //convert into gbk and reformat for BigSCAPE
		//bedops workflow, receive bedfiles with gbk names, and genbank files
		//return gbkfiles in antismash format, return all antismash gbks
		//remove gkb from second bed intersecting, return in separate folder
		//genebank not in antismash, and for this also table with intersection values
		// keep bgc predicted with antismash, rename gbk with 'as' prefix
		// remove bgc predicted with second software intersecting with anitsmash
		// keep names of intersectoin
		//keep gbk predicted with second sofware and not intersecting, save
		// second software folder, rename gbk with 'second software' prefix
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
