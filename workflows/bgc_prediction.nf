/*
General workflow, Take contigs as input and output bgcprediction fasta/gbk
*/

//import modules
include {filterbysize}			from  "../modules/filter"
include {filter_pfamgbk}		from 	"../modules/filter"
include {splitbysize}				from	"../modules/filter"
include {splitgbk}					from	"../modules/bgcgbk"
include {deepbgc_prepare, deepbgc_detect, rungecco, runantismash, runsanntis}	from  "../modules/bgc"
include {prodigal, interproscan, sanntisgbk}  from	"../modules/genepredict"
include {bedops as bedops_sn, bedops as bedops_gc, bedops as bedops_dp}       from  "../modules/catalogue"
include {parse_asresults, parse_snresults, parse_gcresults, parse_dpresults}  from "../modules/processbgc"

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
	sizecontigs=longcontigs.splitFasta(size: "10.MB" ,file:true)
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
		sanntisinput.view()
		runsanntis(sanntisinput)
		gff_sn=runsanntis.out.gff

		parse_snresults(gff_sn, bysizecontigs, gbk)
    bed_sn=parse_snresults.out.bed
    fna_sn=parse_snresults.out.fna
    //reformat sanntis gbk into AS bigscape

    bedops_sn(bed_as.join(bed_sn))
    bed_1=bedops_sn.out.bed
  }
  else{
    bed_1=bed_as
  }

	if(params.gecco){
		//skash/gecco-0.6.3:latest
		rungecco(bysizecontigs) //remove predictions on edges
		gbk_gc=rungecco.out.gbk
		tsv_gc=rungecco.out.tsv

		parse_gcresults(gbk_gc, tsv_gc, bizecontigs)
    bed_gc=parse_gcresults.out.bed
    fna_gc=parsegcresults.out.fna
    //reformat gecco gbk into AS bigscape

    bedops_gc(bed_1.join(bed_gc))
    bed_2=bedops_gc.out.bed
  }
  else{
    bed_2=bed_1
  }

	if(params.deepbgc){
	  deepbgc_prepare(bysizecontigs)
		contig_gbk=deepbgc_prepare.out.gbk
		filter_pfamgbk(contig_gbk)
		gbk=filter_pfamgbk.out.pf_gbk
 		deepbgc_detect(gbk)
		gbk_dp=deepbgc_detect.out.bgc_gbk
    tsv_dp=deepbgc_detect.out.tsv

    parse_dpresults(tsv_dp, gbk_dp)
    bed_dp=parse_dpresults.out.bed
    fna_dp=parse_dp_results.out.fna
    //reformat into Antishmash format

    bedops_dp(bed_2.join(bed_dp))
    bed_3=bedops_dp.out.bed
  }
  else{
    bed_3=bed_2
  }

  //getfasta nucleotides  of BGCs from final beds
  bed_fna=bed_3.join(bysizecontigs)
  //get fna from bed keep sequence name identifier as in bed feature-genebank file
  getbgc_fna(bed_fna)
  bgc_fna=getbgc_fna.out.fna
	//BUILD nucleotide BGCcatalogue,
  //concatenate all BGC_fna into one file
  allbgc_fna=collectFile(name:"allbgc.fna")
  //convert fasta file into mash index
  fna2fnamash(allbgc_fna)
  mashfile=fna2fnamash.out.mash
  //get mash distance between sequences in fna file
  fna_mashtriangle(mashfile)
  fna_distances=fna_mashtirangle.out.list05
  //derreplicate fna all samples, cluster all BGC
  fna_mcl_clust(fna_distances)
  clusters=fna_mcl_clust.out.clusters
  //get representative longest sequence, return list of representative features
  fna_get_representatives(bed_fna, clusters, all_bgcfna)
  bgc_catalogue_bed=fna_get_representatives.out.representative_bed
  bgc_catalogue_fna=fna_get_representatives.out.representative_fna
  //create bowtie index for derreplicated catalogue
  build_index(bgc_catalogue_fna)
  build_genomebed(bgc_catalogue_bed)

  //create bedfile for predicted genes in final bgcs from catalogue (from prodigal or from genebank CDS features)
  //input bgcfna run prodigal and predict genes, output genes
  prodigal_bgc(allbgc_fna)
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
