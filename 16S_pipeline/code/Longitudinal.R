args<-commandArgs(trailingOnly = TRUE)
print(args)

wd<-args[1]
otufile<-args[2]
mths<-args[3]
treefile<-args[4]
dir.create(wd,recursive=T)
setwd(wd)
phenofile<-args[5]
tt_db<-args[6]
input_rank<-args[7]
is_cont<-args[8]
l_size_cutoff<-args[9]
Pheno_arg<-args[10]
Cov_arg<-args[11]
Zero_cov_arg<-args[12]
Time_arg<-args[13]

# print(tt_db)

# wd<-"/home/examples/Longitudinal/result"
# otufile<-"/home/examples/Longitudinal/data/otu_table_mc2_new_ez.txt"
# # mths<-args[3]
# treefile<-"/usr/local/16S_pipeline/file/eztaxon_sina_fastmp.tre"
# dir.create(wd,recursive=T)
# setwd(wd)
# phenofile<-"/home/examples/Longitudinal/data/pheno_togo.csv"
# tt_db<-"ez"
# input_rank<-5
# is_cont<-1
# Pheno_arg<-"BinTrait"
# Cov_arg<-"Cov1,Time"
# Zero_cov_arg<-NULL
# Time_arg<-"Time"
# l_size_cutoff<-2000
# .libPaths("~/tmp")

library(ape)
library(data.table)
library(stringr)
library(nlme)
library(NBZIMM)

# phenofile<-"/home/examples/Longitudinal/data/pheno_togo.csv"
# otufile<-"/home/examples/Longitudinal/data/otu_table_mc2_new_ez.txt"
# treefile<-"/usr/local/16S_pipeline/file/eztaxon_sina_fastmp.tre"
# wd <-"/home/examples/Longitudinal/result"

scale<-1.2

makeCovMat<-function(covariates,data){
  cov_mat<-data[covariates]
  indList_num<-which(sapply(cov_mat,is.numeric))
  indList_factor<-which(!sapply(cov_mat,is.numeric))
  covCateList<-lapply(indList_factor,function(ind_factor){
    fastDummies::dummy_cols(cov_mat[ind_factor])[,-1]
  })
  covCate<-do.call("cbind",covCateList)
  if(is.null(covCate)){
    covResult<-cov_mat[indList_num]
  }else{
    covResult<-cbind(cov_mat[indList_num],covCate)  
  }
  return(covResult)
}

log_print<-function(message="",header="[PROCESS]",method=NULL,data_label=NULL,follwing=NULL){
  Method_DB<-paste(c(method,data_label),collapse=", ")
  if(nchar(Method_DB)==0){
    EnvINfo<-NULL
  }else{
    EnvINfo<-paste0(" (",Method_DB,")\n")
  }
  cat(paste0("\n",header," ",message,"\n",EnvINfo))
  if(!is.null(follwing)) print(follwing)
}

build_formula<-function(Core="Y.tran ~ pheno",Cov=NULL,additional=NULL){
  if(is.null(Cov)){
    addmdl<-""
  }else{
    addmdl<-paste0("+",paste(paste0("Cov.",Cov),collapse="+"))
  }
  if(is.null(additional)){
    mdl<-formula(paste0(Core,addmdl))
  }else{
    mdl<-formula(paste0(Core,addmdl,additional))
  }
}


#### Data load

# path_data<-"http://viva1109.iptime.org/data/"
# phenofile<-file.path(path_data,"temporal_spatial_meta.csv")
# otufile<-file.path(path_data,"temporal_spatial_otuTable.csv")
# taxonomyfile<-file.path(path_data,"temporal_spatial_taxonomy.csv")
# treefile<-file.path(path_data,"temporal_spatial_tree.rds")





# path_data<-"~/analysis/20230201AMAA/"
# phenofile<-file.path(path_data,"data/pheno_togo.csv")
# otufile<-"/data/MICROBIOME/mdprofile/5.ansan/open/ez/otu_table_mc2_new.txt"
# treefile<-"/data/MICROBIOME/mdprofile/5.ansan/open/ez/rep_set.tre"


# phenofile<-"/home/Longitudinal/data/pheno_togo.csv"
# otufile<-"/home/Longitudinal/data/otu_table_mc2_new_ez.txt"
# treefile<-"/usr/local/16S_pipeline/file/eztaxon_sina_fastmp.tre"
# wd <-"/home/Longitudinal/result"
# setwd(wd)

#Check the file paths
log_print("[DESCRIPTION]",message = "Input paths: ",follwing = c(phenofile,otufile,treefile))


##### Metagenome file
log_print(message = "Read data")

data_read_bf<-data.frame(fread(otufile),check.names = F)

data_read<-data_read_bf[,-dim(data_read_bf)[2]]
data_input<-data_read[,-1]
rownames(data_input)<-data_read[,1]
data_input[1:3,1:6]

##### Taxonomy information file

tax_raw<-data_read_bf[,dim(data_read_bf)[2]]
tmp_tax<-data.frame(do.call("rbind",str_split(tax_raw,"; ")),stringsAsFactors = F)
names(tmp_tax)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxonomy_Table<-cbind(X.OTU.ID=data_read_bf[,1],tmp_tax,stringsAsFactors = F)

log_print("[DESCRIPTION]",message = "Heads of taxonomy table: ",follwing = head(taxonomy_Table,3))


##### Phenotype and covariates file
#Phenotype file. Note that SampleID column should be a subset of the second to column names of data_input.
metadata_bf<-data.frame(fread(phenofile),stringsAsFactors=F,check.names = F)

log_print("[DESCRIPTION]",message = "Heads of metadata: ")
print(head(metadata_bf,3))



##### Phylogenetic tree
#Tree file.
# download.file(treefile, "./temporal_spatial_tree.rds")
# tree_raw<-readRDS("./temporal_spatial_tree.rds")
log_print(message = "Read tree file")
tree_raw<-read.tree(treefile)

### Data preparation
#### Basic setting


# Set analysis 

Taxonomy_level<-c("Phylum","Class","Order","Family","Genus","Species")[as.numeric(input_rank)]
Phenotype<-Pheno_arg
ID<-"SubjectID"
Cov<-strsplit(Cov_arg,",")[[1]]
Time<-Time_arg

# OTU filtering - proportion of samples with zero count.
taxon_filter_zero_prop<-0.25 

# Library size
total_readcount<-apply(data_input,2,sum)

#### Match sample

# Match samples to have same order for both of metadata and microbiome data.
ind_sample<-match(metadata_bf[,"SampleID"],colnames(data_input))
metadata_bf$total_readcount<-total_readcount[ind_sample]
otu_df<-data_input[,ind_sample]
otu_mat_bf = as.matrix(otu_df)

passed_ind<-metadata_bf$total_readcount > l_size_cutoff

metadata<-metadata_bf[passed_ind,]
otu_mat<-otu_mat_bf[,passed_ind]

#### OTU filtering
log_print(message = "OTU/ASV filtering")


# calculate non-zero proporiton
non.zero.pr = apply(otu_mat, 1, function(x) {length(x[x != 0])/length(x)} )
df_non.zero.pr = data.frame(zero.p=sort(non.zero.pr, decreasing = T))

# non-zero samples need to be more than 25% of total samples.

log_print("[DESCRIPTION]",message = "OTU/ASV remained: ")

ind_otu_remain<-df_non.zero.pr$zero.p>taxon_filter_zero_prop
sum(ind_otu_remain)
OTU_id_list<-rownames(df_non.zero.pr)[ind_otu_remain]

# OTU Table with remaining samples.
otu_mat_sel = t(otu_mat[OTU_id_list, ])

#### Vector and matrix setting for analysis

y_vec<-as.numeric(metadata[,Phenotype])
id_vec<-as.character(metadata[,ID])
trc_vec<-metadata$total_readcount
time_vec<-metadata[,Time]
covResult<-makeCovMat(Cov,metadata)

#### Setting for the analysis across all taxa

Taxo_togo<-taxonomy_Table[match(OTU_id_list,taxonomy_Table[,1]),]
cnt_taxon<-table(Taxo_togo[,Taxonomy_level])

# Taxa to analize with the order of the number of species taxon it contains.
taxon_anal_list<-names(sort(cnt_taxon,decreasing = T))

# Prune tree with OTUs to analysis.
if(!is.rooted(tree_raw)){
  tree_ori<-root(tree_raw, 1, r = TRUE)
}else{
  tree_ori<-tree_raw
}
tree_chosen<-drop.tip(tree_ori,setdiff(tree_ori$tip.label,OTU_id_list))

# Get counts and subtrees for each taxon.
taxon_analysis_list<-lapply(1:length(taxon_anal_list),function(j){
  target_taxon<-taxon_anal_list[j]
  col_ind<-which(Taxo_togo[,Taxonomy_level]==target_taxon)
  
  if(length(col_ind)==1){
    INDI_SINGLE<-T
  }else{
    INDI_SINGLE<-F
  }
  X_P<-matrix(otu_mat_sel[, col_ind],ncol=length(col_ind))
  colnames(X_P)<-colnames(otu_mat_sel)[col_ind]
  
  # tree_chosen
  if(is.rooted(tree_chosen)){
    tree_chosen_1<-tree_chosen
  }else{
    tree_chosen_1<-multi2di(tree_chosen)
  }
  
  if (length(intersect(colnames(X_P),tree_chosen_1$tip.label))>0 ){
    tree_chosen2_bf<-drop.tip(tree_chosen_1,setdiff(tree_chosen_1$tip.label,colnames(X_P)))
    if(!is.rooted(tree_chosen2_bf)){
      tree_chosen2<-root(tree_chosen2_bf, 1, r = TRUE)
    }else{
      tree_chosen2<-tree_chosen2_bf
    }
  }else{
    tree_chosen2<-NULL
  }

  if(INDI_SINGLE){
    taxon_vec<-X_P
    dim(taxon_vec)<-NULL
  }else{
    taxon_vec<-apply(X_P,1,sum)
  }
  
  Counts_C<-matrix(trc_vec-taxon_vec,ncol=1)
  list(taxon_vec=taxon_vec,tree=tree_chosen2)
})


### Conduct Analysis

#### Run LMM - log transformation


# taxon_vec<-taxon_analysis_list[[1]]$taxon_vec

method<-"LMM_LOG"

log_print(message = "Association analysis",method=method)
# Resonse: log((taxon_vec/trc_vec)+1) 1 is added to avoid log(0)



# random=~1| SID means samples with same SID are correlated with a certain correlation value

resultLMM_LOG<-sapply(taxon_analysis_list,function(input){
  taxon_vec<-input$taxon_vec
  tdata <- data.frame(Y.tran=log((taxon_vec/trc_vec)+1),pheno=y_vec,SID=id_vec,Cov=covResult)
  mdl<-build_formula("Y.tran~pheno",Cov = Cov)
  lme.fit <- try(lme(mdl ,random=~1| SID, data = tdata))
  
  if(class(lme.fit)!="try-error"){
    coef.mat <- summary(lme.fit)$tTable["pheno",]
    LMM_result <- coef.mat
  }else{
    LMM_result<-rep(NA,5)
  }
  return(LMM_result)
})

df_result_LMM_LOG<-cbind(Taxon=taxon_anal_list,data.frame(t(resultLMM_LOG)))
df_result_LMM_LOG$FDR_Corrected_pvalue<-p.adjust(df_result_LMM_LOG$p.value,method="fdr")
write.csv(df_result_LMM_LOG,paste0("Result_",method,".csv"))
write.csv(df_result_LMM_LOG[,c("Taxon","FDR_Corrected_pvalue")],paste0("Pval_",method,".csv"))

log_print(message = "Association analysis done",method=method)

method<-"LMM_ARCSROOT"
log_print(message = "Association analysis",method=method)



resultLMM_asroot<-sapply(taxon_analysis_list,function(input){
  taxon_vec<-input$taxon_vec
  tdata <- data.frame(Y.tran=asin(sqrt(taxon_vec/trc_vec)),pheno=y_vec,SID=id_vec,Cov=covResult)
  mdl<-build_formula("Y.tran~pheno",Cov = Cov)
  
  lme.fit <- try(lme(mdl,random=~1| SID, data = tdata))
  
  if(class(lme.fit)!="try-error"){
    coef.mat <- summary(lme.fit)$tTable["pheno",]
    LMM_result <- coef.mat
  }else{
    LMM_result<-rep(NA,5)
  }
  return(LMM_result)
})

df_result_LMM_asroot<-cbind(Taxon=taxon_anal_list,data.frame(t(resultLMM_asroot)))
df_result_LMM_asroot$FDR_Corrected_pvalue<-p.adjust(df_result_LMM_asroot$p.value,method="fdr")
write.csv(df_result_LMM_asroot,"Result_LMM_arcsineroot.csv")
write.csv(df_result_LMM_asroot[,c("Taxon","FDR_Corrected_pvalue")]," Pval_LMM_arcsineroot.csv")
log_print(message = "Association analysis done",method=method)


method<-"NBMM"

log_print(message = "Association analysis",method=method)



# Note that: offset(log(TRCzinb)), random = ~1|id_vec
# offset if for the adjustment of library size.
mdl<-build_formula("Responsezinb ~ pheno",Cov=Cov,additional = "+ offset(log(TRCzinb))")

NBMMresult<-sapply(taxon_analysis_list,function(input){
  taxon_vec<-input$taxon_vec
  zinbDF<-data.frame(Responsezinb=taxon_vec,pheno=y_vec, IDzinb=id_vec,TRCzinb=trc_vec,Cov=covResult)
  
  NBMM = glmm.nb(mdl, random = ~1|id_vec, verbose = F,data=zinbDF)
  summary(NBMM)$tTable

  if(class(NBMM)[1]!="try-error"){
    NBMMresult= summary(NBMM)$tTable["pheno", ]
  }else{
    NBMMresult<-rep(NA,5)
  }
  return(NBMMresult)
})

df_result_NBMM<-cbind(Taxon=taxon_anal_list,data.frame(t(NBMMresult)))
df_result_NBMM$FDR_Corrected_pvalue<-p.adjust(df_result_NBMM$p.value,method="fdr")
write.csv(df_result_NBMM,"Result_NBMM.csv")
write.csv(df_result_NBMM[,c("Taxon","FDR_Corrected_pvalue")],"Pval_NBMM.csv")
log_print(message = "Association analysis done",method=method)

method<-"ZINBMM"
log_print(message = "Association analysis",method=method)


# Note that: offset(log(TRCzinb)), random = ~1|id_vec
# offset if for the adjustment of library size.
# Nomal case: zero inflation is modeled as same for all the samples.
mdl<-build_formula("Responsezinb ~ pheno",Cov=Cov,additional = "+ offset(log(TRCzinb))")


resultFZINBMM<-sapply(taxon_analysis_list,function(input){
  taxon_vec<-input$taxon_vec
  zinbDF<-data.frame(Responsezinb=taxon_vec,pheno=y_vec, IDzinb=id_vec,TRCzinb=trc_vec,Cov=covResult)
  
  f3 = try(glmm.zinb(mdl, random = ~ 1|IDzinb,data=zinbDF))

  if(class(f3)[1]!="try-error"){
    f4<-summary(f3)
    res_FZINBMM<-f4$tTable[2,]
  }else{
    res_FZINBMM<-rep(NA,5)
  }
  names(res_FZINBMM)<-c("Value", "Std.Error", "DF",     "t-value"   ,   "p-value")
  return(res_FZINBMM)
})

df_resultFZINBMM<-cbind(Taxon=taxon_anal_list,data.frame(t(resultFZINBMM)))
df_resultFZINBMM$FDR_Corrected_pvalue<-p.adjust(df_resultFZINBMM$p.value,method="fdr")
write.csv(df_resultFZINBMM,"Result_ZINBMM.csv")
write.csv(df_resultFZINBMM[,c("Taxon","FDR_Corrected_pvalue")],"Pval_ZINBMM.csv")
log_print(message = "Association analysis done",method=method)

method<-"mTMAT"

log_print("Association analysis",method=method)
source("http://viva1109.duckdns.org/RFunctions_docker/mTMAT.R",encoding = "UTF-8")

output_mTMAT<-mTMAT(Data_Read=otu_df,Meta_Table=metadata,Phenotype=Phenotype,ID=ID,SampleID="SampleID",Time=Time,Cov=Cov,Tree=tree_raw,Taxonomy_Table=taxonomy_Table)

# Significant result

write.csv(output_mTMAT$resultTable,paste0("Result_",method,".csv"))
saveRDS(output_mTMAT,paste0("output_",method,".rds"))
outout_formatted<-output_mTMAT$resultTable[,c(Taxonomy_level,"mTMAT_FDR")]
names(outout_formatted)<-c("Taxon","FDR_Corrected_pvalue")
write.csv(outout_formatted,paste0("Pval_",method,".csv"))

log_print(message = "Association analysis done",method=method)


