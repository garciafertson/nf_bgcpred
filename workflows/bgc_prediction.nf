/*
General workflow, Take contigs as input and output bgcprediction fasta/gbk
*/

//import modules
include {filterbysize}			from  "../modules/filter"
include {filter_pfamgbk}		from 	"../modules/filter"
include {splitbysize}				from	"../modules/filter"
include {splitgbk}					from	"../modules/bgcgbk"
include {deepbgc_prepare; deepbgc_detect; rungecco; runantismash; runsanntis}	from  "../modules/bgc"
include {prodigal; interproscan; sanntisgbk}  from	"../modules/genepredict"
include {bedops_sn}    from    "../modules/catalogue"  
include {bedops_gc}    from    "../modules/catalogue"
include {bedops_dp}    from    "../modules/catalogue" 
include {getbgc_fna; getbgc_bed}	         from    "../modules/catalogue"  
include { mashtriangle; mcl_clust; fna_get_representatives; build_index}	from "../modules/catalogue"
include {build_genomebed; prodigal_bgc; build_bedbgc}	from "../modules/catalogue"
include {parse_asresults; parse_snresults; parse_gcresults; parse_dpresults}  from "../modules/processbgc"

//run metagenomic assembly pipeline using megahit
process modify_contigid_splitfas {
input:
  tuple val(x), path(splitcontig)
output:
  tuple val(baseName),  path(splitcontig), emit: splitcontigs
script:
  split_name=splitcontig.getName()
  baseName = split_name.substring(0, split_name.lastIndexOf('.'))
"""
echo ${baseName}
"""
}



workflow BGCPRED {
	folder=params.contig_folder
	fasta=Channel.fromPath(["${folder}/*contigs.fa", "${folder}/*.fasta", "${folder}/*.fna"])
	filterbysize(fasta)
	longcontigs=filterbysize.out.contigs
	sizecontigs=longcontigs.splitFasta(size: "6.MB" ,file:true)
	modify_contigid_splitfas(sizecontigs)
	bysizecontigs=modify_contigid_splitfas.out.splitcontigs
	// Split large contig file into ~100 MB files
	// Bedfiles chomosome name, start, end, feature name

	runantismash(bysizecontigs)
	gbk_as=runantismash.out.gbk
	parse_asresults(gbk_as)
	bed_as=parse_asresults.out.bed
	fna_as=parse_asresults.out.fna

	if(params.sanntis){
		prodigal(bysizecontigs)
		genesfaa=prodigal.out.genesfaa
		intogbk=bysizecontigs.join(genesfaa)
		sanntisgbk(intogbk)
		gbk=sanntisgbk.out.gbk
		interproscan(genesfaa)
		gff3=interproscan.out.iptsv
		sanntisinput=gff3.join(gbk)
		runsanntis(sanntisinput)
		gff_sn=runsanntis.out.gff
		sn_parse=gff_sn.join(bysizecontigs)
		parse_snresults(sn_parse)
		bed_sn=parse_snresults.out.bed
		fna_sn=parse_snresults.out.fna
 	 }
	else{
		bed_sn=bed_as
	}

	if(params.gecco){
		//skash/gecco-0.6.3:latest
		rungecco(bysizecontigs) //remove predictions on edges
		gbk_gc=rungecco.out.gbk
		tsv_gc=rungecco.out.tsv
		gc_parse=gbk_gc.join(tsv_gc)
		parse_gcresults(gc_parse)
  		bed_gc=parse_gcresults.out.bed
  		fna_gc=parse_gcresults.out.fna
  	}
 	else{
  		bed_gc=bed_as
  	}

	if(params.deepbgc){
	  deepbgc_prepare(bysizecontigs)
		contig_gbk=deepbgc_prepare.out.gbk
		filter_pfamgbk(contig_gbk)
		gbk=filter_pfamgbk.out.pf_gbk
 		deepbgc_detect(gbk)
		gbk_dp=deepbgc_detect.out.bgc_gbk
		tsv_dp=deepbgc_detect.out.tsv
		dp_parse=tsv_dp.join(gbk_dp)
		parse_dpresults(dp_parse)
		bed_dp=parse_dpresults.out.bed
		fna_dp=parse_dpresults.out.fna
		//reformat into Antishmash format
	}
	else{
		bed_dp=bed_as
	}
  bed_snin=bed_as.join(bed_sn)
  bedops_sn(bed_snin)
  bed_snout=bedops_sn.out.bed
  bed_gcin=bed_snout.join(bed_gc)

  bedops_gc(bed_gcin)
  bed_gcout=bedops_gc.out.bed
  bed_dpin=bed_gcout.join(bed_dp)

  bedops_dp(bed_dpin)
  bed_final=bedops_dp.out.bed

  //getfasta nucleotides  of BGCs from final beds
  bed_fna=bed_final.join(bysizecontigs)
  //get fna from bed keep sequence name identifier as in bed feature-genebank file
  getbgc_fna(bed_fna)
  bgc_fna=getbgc_fna.out.fna
  //BUILD nucleotide BGCcatalogue,
  //concatenate all BGC_fna into one file
  bgc_fna.view()
  allbgc_fna=bgc_fna.collectFile(name:"allbgc.fna")
  
  //get mash distance between sequences in fna file
  mashtriangle(allbgc_fna)
  fna_distances=mashtriangle.out.list05
  
  //derreplicate fna all samples, cluster all BGC
  mcl_clust(fna_distances)
  clusters=mcl_clust.out.clusters
  //get representative longest sequence, return list of representative features

  getbgc_bed(bed_final)
  allbgc_bed=getbgc_bed.out.bed.collectFile(name:"allbgc.bed")

  fna_get_representatives(allbgc_bed, clusters, allbgc_fna)
  bgc_catalogue_bed=fna_get_representatives.out.representative_bed
  bgc_catalogue_fna=fna_get_representatives.out.representative_fna
  //create bowtie index for derreplicated catalogue
  build_index(bgc_catalogue_fna)
  build_genomebed(bgc_catalogue_bed)
  
  //create bedfile for predicted genes in final bgcs from catalogue (from prodigal or from genebank CDS features)
  //input bgcfna run prodigal and predict genes, output genes
  prodigal_bgc(bgc_catalogue_fna)
  allbgc_genes=prodigal_bgc.out.gff
  //input prodigal gff(gbk) and reformat into bed file.
  build_bedbgc(allbgc_genes)

  ////////////////////
  //recover genebank files from representatives for bigscape
  //antismash genebank ready for bigscape,
  //gecco genebank, convert into AS format for bigscape
  //deepbgc genbank, convert into AS format for bigscape
  //sanntis, convert fna into genebank using prodigal and convert into AS format

	//*** compare final predicted sequecnes and keep:
	// Check 1 - all AS predictions
	// * missing 2 - only intersections with other predictiors >1
	// Check 3 - only very high scoring predicions single predictor
}
