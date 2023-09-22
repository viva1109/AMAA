library(dada2)
library(phyloseq)

###########################################################################
silvaTsetFile = "/data/sharedcode/scpark/bash/16S_pipeline/file/silva_nr_v132_train_set.fa.gz"
silvaAssFile = "/data/sharedcode/scpark/bash/16S_pipeline/file/silva_species_assignment_v132.fa.gz"

ezTsetFile = "/data/sharedcode/scpark/bash/16S_pipeline/file/mod_eztaxon_qiime_full.fasta.gz"
ezAssFile  = "/data/sharedcode/scpark/bash/16S_pipeline/file/mod_eztaxon_qiime_full.fasta.species.gz"

args = commandArgs(trailingOnly = TRUE)
trimdir = args[1]
forid   = args[2]
revid   = args[3]
dbid    = args[4]
ncore   = as.numeric(args[5])
sstep   = as.numeric(args[6])
alloutdir  = args[7]
fMaxEE   = as.numeric(args[8])
rMaxEE   = as.numeric(args[9])
fTL      = as.numeric(args[10])
rTL      = as.numeric(args[11])

print(paste0("dada2 parameters"))
print(paste0("  trimdir: ", trimdir))
print(paste0("  forid  : ", forid))
print(paste0("  revid  : ", revid))
print(paste0("  dbid   : ", dbid))
print(paste0("  ncore  : ", ncore))
print(paste0("  step   : ", sstep))
print(paste0("  outdir : ", alloutdir))
print(paste0("  maxEE(for)    : ", fMaxEE))
print(paste0("  maxEE(rev)    : ", rMaxEE))
print(paste0("  truncLen(for) : ", fTL))
print(paste0("  truncLen(rev) : ", rTL))

if(length(args) != 7) quit()

# trimPath = "/data/sharedcode/scpark/bash/16S_pipeline/testdir/1.trimmed"
# forid = "_R1_001"
# revid = "_R2_001"

executedada = function(trimdir,outdir,fMaxEE=3,rMaxEE=5,fTL=275,rTL=250,ncore){
    fnFs = sort(list.files(trimdir, pattern = paste0(forid,".fastq"),  full.names = TRUE))
    fnRs = sort(list.files(trimdir, pattern = paste0(revid,".fastq") , full.names = TRUE))
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    print(paste0("Total pairs : ",length(fnFs)))
    if(length(fnFs) == 0) {
        print("No file in trimmed directory")
        quit()
    }
    dfiltpath = file.path(dirname(trimdir), "dada_filt")
    if(!file_test("-d", dfiltpath)) dir.create(dfiltpath)
    filtFs = file.path(dfiltpath, paste0(sample.names, paste0("_filt",forid,".fastq")))
    filtRs = file.path(dfiltpath, paste0(sample.names, paste0("_filt",revid,".fastq")))
    print(paste("Total filtered file : ",length(filtFs),length(filtRs)))
    filtRes = filterAndTrim(fnFs, filtFs,fnRs, filtRs, maxEE=c(fMaxEE,rMaxEE),compress=F,
                            truncLen=c(fTL,rTL),multithread=ncore, verbose=TRUE)

    set.seed(100)
    errF <- learnErrors(filtFs, nbases=1e8, multithread=ncore)
    errR <- learnErrors(filtRs, nbases=1e8, multithread=ncore)
    
    mergers <- vector("list", length(sample.names))
    names(mergers) <- sample.names
    
    idx = 1
    for(sam in sample.names) {
        cat("Processing:", sam, "\n")
        derepF = derepFastq(filtFs[idx])
        ddF = dada(derepF, err=errF, multithread=ncore)
        derepR = derepFastq(filtRs[idx])
        ddR = dada(derepR, err=errR, multithread=ncore)
        merger = mergePairs(ddF, derepF, ddR, derepR)
        mergers[[sam]] <- merger
        idx = idx +1
    }
    
    seqtab = makeSequenceTable(mergers)
   
    if(!file_test("-d", outdir)) dir.create(outdir)
    seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=ncore)
    rdsFile = file.path(outdir,"seqtab.nochim.rds")
    saveRDS(seqtab.nochim, rdsFile)
    print(paste0("dada step 1 done"))
}
    
executeTaxAssign = function(rdsFile,db,outdir,ncore){
    if(!db %in% c("gg","ez","silva")){
        print(paste0("No such db : ",db))
        quit()
    }
    if(!file.exists(rdsFile)){
        print(paste0("No seqchim rds file , rerun from start"))
        quit()
    }
    seqtab.nochim = readRDS(rdsFile)
    # Assign taxonomy
    if( db == "silva") {
        trainSetFile = silvaTsetFile
        assignFile = silvaAssFile
    }else if(db == "ez") {
        trainSetFile = ezTsetFile
        assignFile = ezAssFile
    }
    tax = assignTaxonomy(seqtab.nochim, trainSetFile, multithread=ncore)
    if(db != "ez") {
        tax = addSpecies(tax, assignFile)
    }else{
        tax = addSpecies(tax, assignFile,allowMultiple = 1)
    }
    
    # Write to disk
    if(!dir.exists(outdir)) dir.create(outdir)
    saveRDS(seqtab.nochim,  file.path(outdir,paste0(db,"_seqtab_nochim.rds")))
    saveRDS(tax, file.path(outdir,paste0(db,"_tax_final.rds")))
    
    ps = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(tax))
    dna = Biostrings::DNAStringSet(taxa_names(ps))
    names(dna) = taxa_names(ps)
    ps = merge_phyloseq(ps, dna)
    taxa_names(ps) = paste0("ASV", seq(ntaxa(ps)))
    
    OTU = as(otu_table(ps), "matrix")
    if(!taxa_are_rows(ps)){OTU <- t(OTU)}
    TAX = as(tax_table(ps),"matrix")
    full_tax =  apply( TAX[ , colnames(TAX) ] , 1 , paste0 , collapse = ";" )
    
    fOTU = cbind(OTU,taxonomy=full_tax)
    otuFile = file.path(outdir,"otu_table_w.txt")
    write.table(fOTU,otuFile,sep = "\t",quote = F,col.names = NA)
    print(paste0("dada step 2 done"))
    return(otuFile)
}

trimPath = normalizePath(trimdir)
if(sstep<2){
    executedada(trimdir=trimPath,outdir = alloutdir,
                fMaxEE=fMaxEE,rMaxEE =rMaxEE,fTL = fTL, rTL = rTL,ncore = ncore)
}
rdsFile = file.path(alloutdir,"seqtab.nochim.rds")
executeTaxAssign(rdsFile,dbid,alloutdir,ncore)

