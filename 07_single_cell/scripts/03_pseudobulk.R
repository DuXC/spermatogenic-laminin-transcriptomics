#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE,width=160)
.libPaths(c(normalizePath('env/R_library'),.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
set.seed(20260909)
S <- '07_single_cell'
write_tsv <- function(x,path,row.names=FALSE) {
  if(grepl('\\.gz$',path)) {con<-gzfile(path,'wt');on.exit(close(con))} else con<-path
  write.table(x,con,sep='\t',quote=FALSE,row.names=row.names,na='NA')
}
sets <- fromJSON('02_annotation/fixed_gene_sets.json')$primary
samples <- read.delim(file.path(S,'03_processed/pseudobulk_samples.tsv'),check.names=FALSE)
counts <- as.matrix(read.delim(file.path(S,'03_processed/pseudobulk_counts.tsv.gz'),row.names=1,check.names=FALSE))
stopifnot(identical(colnames(counts),samples$sample_id),all(counts>=0),all(counts==floor(counts)))
stats <- list(); pathways<-list(); qc<-list(); coverage<-list()
for (variant in unique(samples$variant)) {
 for (typ in unique(samples$cell_type)) {
  s<-samples[samples$variant==variant & samples$cell_type==typ & samples$technology=='BD_Rhapsody',]
  ok<-nrow(s)==6 && all(s$n_cells>=30) && sum(s$group=='OA')==3 && sum(s$group=='iNOA')==3 && all(s$annotation_confidence=='high')
  coverage[[length(coverage)+1]]<-data.frame(variant=variant,cell_type=typ,n_OA=sum(s$group=='OA' & s$n_cells>=30),n_iNOA=sum(s$group=='iNOA' & s$n_cells>=30),eligible=ok)
  if(!ok) next
  s<-s[order(factor(s$group,levels=c('OA','iNOA')),s$donor_id),]
  x<-counts[,s$sample_id,drop=FALSE]
  y<-DGEList(x)
  keep<-rowSums(cpm(y)>=1)>=3
  y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
  s$group<-factor(s$group,levels=c('OA','iNOA'))
  design<-model.matrix(~group,s)
  v<-voom(y,design,plot=FALSE,save.plot=TRUE)
  fit<-eBayes(lmFit(v,design),robust=TRUE)
  tab<-topTable(fit,coef=2,number=Inf,sort.by='none',confint=TRUE)
  tab$gene<-rownames(tab);tab$cell_type<-typ;tab$variant<-variant
  tab$contrast<-'iNOA_minus_OA';tab$n_OA<-3;tab$n_iNOA<-3
  stats[[length(stats)+1]]<-tab
  idx<-ids2indices(sets,rownames(y))
  cam<-camera(v,idx,design,contrast=2,inter.gene.cor=NA,allow.neg.cor=FALSE,sort=FALSE)
  cam$pathway<-rownames(cam);cam$cell_type<-typ;cam$variant<-variant
  cam$n_OA<-3;cam$n_iNOA<-3
  # Descriptive average of all measured member-gene log2 effects; camera supplies competitive inference.
  cam$mean_member_logFC<-vapply(idx,function(j)mean(fit$coefficients[j,2]),numeric(1))
  pathways[[length(pathways)+1]]<-cam
  lib<-y$samples;lib$sample_id<-rownames(lib);lib$cell_type<-typ;lib$variant<-variant;lib$n_genes<-nrow(y)
  qc[[length(qc)+1]]<-lib
  fn<-paste0(variant,'_',typ)
  write_tsv(data.frame(gene=rownames(v$E),v$E,check.names=FALSE),file.path(S,'04_results',paste0('voom_logCPM_',fn,'.tsv.gz')))
  png(file.path(S,'04_results',paste0('voom_mean_variance_',fn,'.png')),width=1000,height=700,res=130)
  plot(v$voom.xy$x,v$voom.xy$y,pch='.',xlab='Average log-count size + 0.5',ylab='Sqrt residual standard deviation',main=paste(typ,variant))
  lines(v$voom.line,col='#B45643',lwd=2);dev.off()
  cat('DONE',variant,typ,nrow(y),'genes\n')
 }
}
write_tsv(do.call(rbind,coverage),file.path(S,'04_results/pseudobulk_eligibility.tsv'))
if(length(stats)) {
  tab<-do.call(rbind,stats);rownames(tab)<-NULL
  write_tsv(tab,file.path(S,'04_results/pseudobulk_all_genes.tsv.gz'))
  cam<-do.call(rbind,pathways);rownames(cam)<-NULL
  cam$FDR_all_types_six_sets<-ave(cam$PValue,cam$variant,FUN=function(p)p.adjust(p,'BH'))
  write_tsv(cam,file.path(S,'04_results/pseudobulk_camera.tsv'))
  write_tsv(do.call(rbind,qc),file.path(S,'04_results/pseudobulk_library_QC.tsv'))
}
sink(file.path(S,'00_admin/pseudobulk_R_sessionInfo.txt'));sessionInfo();sink()
