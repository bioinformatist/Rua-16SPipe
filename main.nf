#!/usr/bin/env nextflow

// ANSI_RESET = "\u001B[0m"
// ANSI_BLACK = "\u001B[30m"
// ANSI_RED = "\u001B[31m"
// ANSI_GREEN = "\u001B[32m"
// ANSI_YELLOW = "\u001B[33m"
// ANSI_BLUE = "\u001B[34m"
// ANSI_PURPLE = "\u001B[35m"
// ANSI_CYAN = "\u001B[36m"
// ANSI_WHITE = "\u001B[37m"

// def helpMessage() {

//     log.info """
    
//     $ANSI_YELLOW
//     Usage:

//     The typical command for running the pipeline is as follows:
//     $ANSI_RESET
//     nextflow run TODO

//     Mandatory arguments:
//         --reads_mode
    
//     """
//     .stripIndent()
// }

// /*

//  * SET UP CONFIGURATION VARIABLES

//  */


// // Show help emssage
// params.help = false
// if (params.help){
//     helpMessage()
//     exit 0
// }

// params.inputdir = "$baseDir/data/mixed"
params.silva_align = "$baseDir/dbs/silva.nr_v132.align"
params.silva_tax = "$baseDir/dbs/silva.nr_v132.tax"
params.silva_subset = "$baseDir/dbs/silva.v45.fasta"
params.outdir = 'results'
params.mothur_prefix = 'hehe'

// log.info """\

//          Rua-16S    P I P E L I N E    

//          =============================

//          mode         : ${params.reads_mode}   

//          Data         : TODO

//          Database     : ${params.silva_align}

//          Sub-database : ${params.silva_subset}

//          Tax          : ${params.silva_tax}

//          Outdir       : ${params.outdir}

//          """

//          .stripIndent()


if (params.reads_mode == 'mixed') {
    // Channel
    //     .fromFilePairs(params.ffastq, params.rfastq, params.oligo)
    //     .into {ffastq_ch, rfastq_ch, oligo_ch}

    pre_process = """
    make.contigs(ffastq=$launchDir/$params.ffastq, rfastq=$launchDir/$params.rfastq, oligos=$launchDir/$params.oligo, pdiffs=6, bdiffs=1, checkorient=t, processors=56)
    set.dir(input=${file(params.ffastq).getParent()}, output=mothur_temp)
    """
    params.mothur_prefix = file(params.ffastq).getBaseName()
    
}

if (params.reads_mode == 'demuxed') {
    pre_process = """
    make.file(inputdir=$launchDir/$params.inputdir, type=fastq, prefix=$params.mothur_prefix)
    make.contigs(file=${params.mothur_prefix}.files, checkorient=t, processors=50)
    set.dir(input=$launchDir/$params.inputdir, output=mothur_temp)
    """
}

process oneNoodlesDragon {

    publishDir params.outdir, mode: 'copy', saveAs: {"rep.fasta"}

    output:

    file "final.opti_mcc.shared" into shared_file
    file 'final.taxonomy' into tax_file

    file "mothur_temp/${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta" into 'rep.fasta'

    """
    #!/usr/bin/env mothur
    $pre_process
    summary.seqs(fasta=${params.mothur_prefix}.trim.contigs.fasta)

    screen.seqs(fasta=${params.mothur_prefix}.trim.contigs.fasta, group=${params.mothur_prefix}.contigs.groups, summary=${params.mothur_prefix}.trim.contigs.summary, optimize=minlength-maxlength, maxambig=0, criteria=97.5)
    summary.seqs()
    unique.seqs(fasta=${params.mothur_prefix}.trim.contigs.good.fasta)
    count.seqs(name=${params.mothur_prefix}.trim.contigs.good.names, group=${params.mothur_prefix}.contigs.good.groups)
    summary.seqs(count=${params.mothur_prefix}.trim.contigs.good.count_table)
    align.seqs(fasta=${params.mothur_prefix}.trim.contigs.good.unique.fasta, reference=$params.silva_subset)
    summary.seqs(fasta=${params.mothur_prefix}.trim.contigs.good.unique.align, count=${params.mothur_prefix}.trim.contigs.good.count_table)

    screen.seqs(fasta=${params.mothur_prefix}.trim.contigs.good.unique.align, count=${params.mothur_prefix}.trim.contigs.good.count_table, summary=${params.mothur_prefix}.trim.contigs.good.unique.summary, optimize=start-end, criteria=97.5, maxhomop=8, processors=50)
    summary.seqs(fasta=current, count=current)
    filter.seqs(fasta=${params.mothur_prefix}.trim.contigs.good.unique.good.align, vertical=T, trump=.)
    unique.seqs(fasta=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.fasta, count=${params.mothur_prefix}.trim.contigs.good.good.count_table)
    pre.cluster(fasta=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.fasta, count=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.count_table, diffs=3)
    chimera.vsearch(fasta=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
    remove.seqs(fasta=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
    summary.seqs(fasta=current, count=current)
    classify.seqs(fasta=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=$params.silva_align, taxonomy=$params.silva_tax, cutoff=80, iters=1000)
    remove.lineage(fasta=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Eukaryota)
    summary.tax(taxonomy=current, count=current)
        
    cluster.split(fasta=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.03, processors=50)
    make.shared(list=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
    classify.otu(list=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, count=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, label=0.03)

    get.oturep(column=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, name=${params.mothur_prefix}.trim.contigs.good.names, list=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.list, fasta=${params.mothur_prefix}.trim.contigs.good.unique.fasta)

    rename.file(taxonomy=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy, shared=${params.mothur_prefix}.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared, prefix=final)
    """
}

process makeBiom {

    input:

    file shared_file
    file tax_file

    output:

    file 'final.opti_mcc.0.03.biom' into biom_file

    """
    mothur "#make.biom(shared=$shared_file, constaxonomy=$tax_file)"
    """
}

process biom2consensusLineage {

    publishDir params.outdir, mode: 'copy'
    
    input:

    file biom_file

    output:

    file "consensusLineage.txt"

    """
    biom convert -i $biom_file -o consensusLineage.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage" --table-type="OTU table"
    """
}


process rPart {

    publishDir params.outdir, mode: 'copy'

    input:

    file shared_file
    file tax_file

    output:

    file "Genus.RA.xlsx"

    """
    #!/usr/bin/env Rscript
    library(phyloseq)
    library(xlsx)
    library(dplyr)
    mothur.data <- import_mothur(mothur_shared_file = '$shared_file', mothur_constaxonomy_file = '$tax_file')
    colnames(tax_table(mothur.data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    RA <- mothur.data %>%  transform_sample_counts(function(x) {x/sum(x)} ) %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE)
    RA.tax <- cbind(as.data.frame(tax_table(RA)), as.data.frame(otu_table(RA)))
    write.xlsx(RA.tax, 'Genus.RA.xlsx', sheetName = "Genus.RA",
           row.names = TRUE, append = FALSE)
    """
}