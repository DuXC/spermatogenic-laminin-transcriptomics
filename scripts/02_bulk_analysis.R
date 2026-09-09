#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE, width=160)
.libPaths(c(normalizePath('env/R_library'), .libPaths()))
suppressPackageStartupMessages({library(limma);library(jsonlite)})
set.seed(20260909)
stopifnot(dir.exists('04_results'))
write_tsv <- function(x,path,row.names=FALSE) {
  if(is.character(path) && grepl('\\.gz$',path)) {
    con<-gzfile(path,'wt');on.exit(close(con));write.table(x,con,sep='\t',quote=FALSE,row.names=row.names,na='NA')
  } else write.table(x,path,sep='\t',quote=FALSE,row.names=row.names,na='NA')
}
read_tsv <- function(path) read.delim(path,check.names=FALSE)
sets <- fromJSON('02_annotation/fixed_gene_sets.json')
allsets <- c(sets$primary,sets$exploratory_markers)
primary_names <- names(sets$primary)
for (p in c('04_results/primary','04_results/sensitivity','04_results/qc')) dir.create(p,showWarnings=FALSE)

load_cohort <- function(acc,variant='primary') {
  s<-read_tsv(paste0('02_annotation/',acc,'_samples.tsv'))
  a<-read_tsv(paste0('02_annotation/',acc,'_probe_to_current_gene.tsv.gz'))
  x<-as.matrix(read.delim(paste0('03_processed/',acc,'_uploaded_probe_values.tsv.gz'),row.names=1,check.names=FALSE))
  stopifnot(identical(colnames(x),s$gsm),identical(rownames(x),as.character(a$probe_id)))
  keep<-as.character(a$analysis_eligible)=='True'
  stopifnot(sum(keep)>1000)
  x<-x[keep,,drop=FALSE];a<-a[keep,,drop=FALSE]
  if (acc=='GSE4797') {
    floor<-if (variant=='floor01') .1 else .01
    x<-log2(pmax(x,floor))
    if (variant!='uploaded_log_only') x<-normalizeBetweenArrays(x,method='quantile')
    s$group<-factor(s$group,levels=c('JS10','JS8','JS5','JS2'))
  } else s$group<-factor(s$group,levels=c('normal','impaired'))
  # Equal-weight mean of eligible probe log values; group labels do not select probes.
  e<-avereps(x,ID=a$gene)
  stopifnot(!anyNA(e),!anyDuplicated(rownames(e)))
  design<-model.matrix(~0+group,s);colnames(design)<-levels(s$group)
  if(acc=='GSE4797') {
    cm<-makeContrasts(JS8_vs_JS10=JS8-JS10,JS5_vs_JS10=JS5-JS10,JS2_vs_JS10=JS2-JS10,JS5_vs_JS8=JS5-JS8,JS2_vs_JS5=JS2-JS5,levels=design)
  } else cm<-makeContrasts(impaired_vs_normal=impaired-normal,levels=design)
  list(acc=acc,s=s,e=e,design=design,cm=cm,probe=x,annotation=a)
}

limma_tables <- function(d) {
  fit<-eBayes(contrasts.fit(lmFit(d$e,d$design),d$cm),robust=TRUE,trend=TRUE)
  out<-lapply(seq_len(ncol(d$cm)),function(j){
    tab<-topTable(fit,coef=j,number=Inf,sort.by='none',confint=TRUE)
    tab$gene<-rownames(tab);tab$contrast<-colnames(d$cm)[j];tab$cohort<-d$acc
    tab$SE_moderated<-fit$stdev.unscaled[,j]*sqrt(fit$s2.post)
    tab[,c('cohort','contrast','gene','logFC','CI.L','CI.R','SE_moderated','AveExpr','t','P.Value','adj.P.Val','B')]
  })
  out<-do.call(rbind,out);rownames(out)<-NULL
  out$FDR_primary_family<-NA_real_
  idx<-if(d$acc=='GSE4797')out$contrast %in% colnames(d$cm)[1:3] else rep(TRUE,nrow(out))
  out$FDR_primary_family[idx]<-p.adjust(out$P.Value[idx],'BH')
  if(d$acc=='GSE4797') {
    ft<-topTable(fit,coef=1:3,number=Inf,sort.by='none');ft$gene<-rownames(ft)
    write_tsv(ft,'04_results/primary/GSE4797_four_group_F.tsv.gz')
    dt<-cbind(intercept=1,johnsen_decrease=10-d$s$johnsen)
    tr<-eBayes(lmFit(d$e,dt),robust=TRUE,trend=TRUE)
    tt<-topTable(tr,coef=2,number=Inf,sort.by='none',confint=TRUE);tt$gene<-rownames(tt)
    write_tsv(tt,'04_results/primary/GSE4797_johnsen_linear_association.tsv.gz')
  }
  out
}

scores_for <- function(d) {
  z<-t(scale(t(d$e)))
  vals<-lapply(allsets,function(gs) {
    g<-intersect(gs,rownames(z))
    if(length(g)<3) rep(NA_real_,ncol(z)) else colMeans(z[g,,drop=FALSE])
  })
  sc<-do.call(rbind,vals);colnames(sc)<-d$s$gsm
  sc
}
score_effects <- function(d,sc,adjust=FALSE) {
  design<-d$design;cm<-d$cm
  if(adjust) {design<-cbind(design,postmeiotic=as.numeric(sc['Postmeiotic_markers',]));cm<-rbind(cm,postmeiotic=0)}
  inv<-solve(crossprod(design)); df<-nrow(design)-ncol(design)
  out<-list()
  for (nm in rownames(sc)) {
    if(anyNA(sc[nm,]) || (adjust && nm=='Postmeiotic_markers'))next
    f<-lm.fit(design,sc[nm,]);v<-sum(f$residuals^2)/df
    for(j in seq_len(ncol(cm))) {
      est<-sum(cm[,j]*f$coefficients);se<-sqrt(drop(t(cm[,j])%*%inv%*%cm[,j])*v)
      out[[length(out)+1]]<-data.frame(cohort=d$acc,set=nm,contrast=colnames(cm)[j],effect=est,SE=se,CI_low=est-qt(.975,df)*se,CI_high=est+qt(.975,df)*se,P=2*pt(-abs(est/se),df),df=df,adjusted_for_postmeiotic=adjust,design_condition_number=kappa(design))
    }
  }
  o<-do.call(rbind,out);o$FDR_family<-NA_real_
  main<-if(d$acc=='GSE4797')colnames(d$cm)[1:3] else colnames(d$cm)
  for(fam in c('primary','markers','secondary')) {
    idx<-switch(fam,primary=o$set %in% primary_names & o$contrast %in% main,markers=!o$set %in% primary_names & o$contrast %in% main,secondary=!o$contrast %in% main)
    o$FDR_family[idx]<-p.adjust(o$P[idx],'BH')
  }
  o
}
camera_tests <- function(d,cor=.01) {
  ind<-ids2indices(sets$primary,rownames(d$e));ind<-ind[lengths(ind)>=10]
  out<-lapply(seq_len(ncol(d$cm)),function(j) {
    tab<-camera(d$e,index=ind,design=d$design,contrast=d$cm[,j],inter.gene.cor=cor,trend.var=TRUE)
    tab$set<-rownames(tab);tab$contrast<-colnames(d$cm)[j];tab$cohort<-d$acc;tab$cor_mode<-if(is.na(cor))'estimated' else 'fixed_0.01';tab
  })
  o<-do.call(rbind,out);o$FDR_family<-NA_real_
  main<-if(d$acc=='GSE4797')o$contrast %in% colnames(d$cm)[1:3] else rep(TRUE,nrow(o))
  o$FDR_family[main]<-p.adjust(o$PValue[main],'BH')
  o$FDR_family[!main]<-p.adjust(o$PValue[!main],'BH')
  rownames(o)<-NULL;o
}
qc <- function(d) {
  pca<-prcomp(t(d$e),center=TRUE,scale.=FALSE)
  ve<-pca$sdev^2/sum(pca$sdev^2)
  pcs<-cbind(d$s,data.frame(pca$x[,1:5]));write_tsv(pcs,paste0('04_results/qc/',d$acc,'_PCA.tsv'))
  write_tsv(data.frame(PC=seq_along(ve),variance_fraction=ve),paste0('04_results/qc/',d$acc,'_PCA_variance.tsv'))
  r<-cor(d$e,method='pearson');write_tsv(data.frame(gsm=rownames(r),r,check.names=FALSE),paste0('04_results/qc/',d$acc,'_sample_correlations.tsv'))
  sm<-data.frame(gsm=d$s$gsm,group=d$s$group,median_log_value=apply(d$e,2,median),IQR=apply(d$e,2,IQR),median_correlation=sapply(seq_len(nrow(r)),function(i)median(r[i,-i])),within_group_median_correlation=sapply(seq_len(nrow(r)),function(i)median(r[i,which(d$s$group==d$s$group[i] & seq_len(nrow(r))!=i)])))
  write_tsv(sm,paste0('04_results/qc/',d$acc,'_sample_QC.tsv'))
}

datasets<-list();stats<-list();scorelist<-list();coverage<-list();summary<-list()
for(acc in c('GSE4797','GSE145467')) {
  cat('ANALYZE',acc,'\n');d<-load_cohort(acc);datasets[[acc]]<-d
  tab<-limma_tables(d);stats[[acc]]<-tab
  write_tsv(tab,paste0('04_results/primary/',acc,'_limma_all_contrasts.tsv.gz'))
  write_tsv(data.frame(gene=rownames(d$e),d$e,check.names=FALSE),gzfile(paste0('03_processed/',acc,'_gene_log_expression.tsv.gz')))
  qc(d)
  for(nm in names(allsets))coverage[[length(coverage)+1]]<-data.frame(cohort=acc,set=nm,total_genes=length(allsets[[nm]]),measured_genes=length(intersect(allsets[[nm]],rownames(d$e))),coverage_fraction=length(intersect(allsets[[nm]],rownames(d$e)))/length(allsets[[nm]]),genes=paste(intersect(allsets[[nm]],rownames(d$e)),collapse='|'),missing_genes=paste(setdiff(allsets[[nm]],rownames(d$e)),collapse='|'))
  sc<-scores_for(d);scorelist[[acc]]<-sc
  write_tsv(data.frame(set=rownames(sc),sc,check.names=FALSE),paste0('04_results/primary/',acc,'_module_scores.tsv'))
  write_tsv(score_effects(d,sc),paste0('04_results/primary/',acc,'_module_effects.tsv'))
  write_tsv(score_effects(d,sc,adjust=TRUE),paste0('04_results/sensitivity/',acc,'_postmeiotic_adjusted_module_effects.tsv'))
  ca<-camera_tests(d);write_tsv(ca,paste0('04_results/primary/',acc,'_camera.tsv'))
  write_tsv(camera_tests(d,cor=NA_real_),paste0('04_results/sensitivity/',acc,'_camera_estimated_correlation.tsv'))
  for(cn in unique(tab$contrast)){
    tt<-tab[tab$contrast==cn,];summary[[length(summary)+1]]<-data.frame(cohort=acc,contrast=cn,genes=nrow(tt),FDR05=sum(tt$adj.P.Val<.05),up_FDR05=sum(tt$adj.P.Val<.05 & tt$logFC>0),down_FDR05=sum(tt$adj.P.Val<.05 & tt$logFC<0),FDR05_abslogFC1=sum(tt$adj.P.Val<.05 & abs(tt$logFC)>=1),FDR05_primary_family=sum(tt$FDR_primary_family<.05,na.rm=TRUE))
  }
  print(ca[,c('set','contrast','Direction','PValue','FDR_family')])
}
write_tsv(do.call(rbind,coverage),'02_annotation/pathway_coverage.tsv')
write_tsv(do.call(rbind,summary),'04_results/primary/differential_expression_summary.tsv')

# Fixed normalization sensitivity; the retained probe universe stays unchanged.
for(variant in c('floor01','uploaded_log_only')) {
  d<-load_cohort('GSE4797',variant)
  write_tsv(camera_tests(d),paste0('04_results/sensitivity/GSE4797_',variant,'_camera.tsv'))
  write_tsv(score_effects(d,scores_for(d)),paste0('04_results/sensitivity/GSE4797_',variant,'_module_effects.tsv'))
  f<-eBayes(contrasts.fit(lmFit(d$e,d$design),d$cm),robust=TRUE,trend=TRUE)
  base<-stats[['GSE4797']];o<-list()
  for(j in 1:3){t<-topTable(f,coef=j,n=Inf,sort.by='none');b<-base[base$contrast==colnames(d$cm)[j],];b<-b[match(rownames(t),b$gene),];o[[j]]<-data.frame(variant=variant,contrast=colnames(d$cm)[j],spearman_logFC=cor(t$logFC,b$logFC,method='spearman'),same_sign_fraction=mean(sign(t$logFC)==sign(b$logFC)),genes=nrow(t))}
  write_tsv(do.call(rbind,o),paste0('04_results/sensitivity/GSE4797_',variant,'_gene_agreement.tsv'))
}

# Sample influence without data-dependent removal: refit every leave-one-out module model.
loo<-list()
for(acc in names(datasets)) {
  d<-datasets[[acc]];sc<-scorelist[[acc]]
  for(i in seq_len(nrow(d$s))){dd<-d;dd$s<-d$s[-i,];dd$design<-d$design[-i,,drop=FALSE];dd$e<-d$e[,-i,drop=FALSE];o<-score_effects(dd,sc[,-i,drop=FALSE]);o$omitted_gsm<-d$s$gsm[i];loo[[length(loo)+1]]<-o}
}
write_tsv(do.call(rbind,loo),'04_results/sensitivity/leave_one_out_module_effects.tsv')

ecm<-unique(unlist(sets$primary));cons<-list();merged<-list()
b<-stats[['GSE145467']]
for(cn in colnames(datasets[['GSE4797']]$cm)[1:3]){
  a<-stats[['GSE4797']];a<-a[a$contrast==cn,]
  m<-merge(a,b,by='gene',suffixes=c('_GSE4797','_GSE145467'));m$in_fixed_ECM_union<-m$gene %in% ecm;merged[[cn]]<-m
  for(scope in c('all_common_genes','fixed_ECM_union')){
    x<-if(scope=='all_common_genes')m else m[m$in_fixed_ECM_union,]
    cons[[length(cons)+1]]<-data.frame(contrast=cn,scope=scope,genes=nrow(x),same_sign_n=sum(sign(x$logFC_GSE4797)==sign(x$logFC_GSE145467)),same_sign_fraction=mean(sign(x$logFC_GSE4797)==sign(x$logFC_GSE145467)),spearman_logFC=cor(x$logFC_GSE4797,x$logFC_GSE145467,method='spearman'),both_FDR05_same_sign=sum(x$adj.P.Val_GSE4797<.05 & x$adj.P.Val_GSE145467<.05 & sign(x$logFC_GSE4797)==sign(x$logFC_GSE145467)))
  }
}
write_tsv(do.call(rbind,cons),'04_results/primary/cross_cohort_direction_summary.tsv')
write_tsv(do.call(rbind,merged),gzfile('04_results/primary/cross_cohort_gene_effects.tsv.gz'))
# Post-run QA extension: restrict every score to exactly the same measurable genes.
# This checks platform coverage and does not replace the frozen primary analyses.
original_allsets<-allsets
common<-intersect(rownames(datasets[['GSE4797']]$e),rownames(datasets[['GSE145467']]$e))
allsets<-lapply(allsets,intersect,y=common)
shared_scores<-list()
for(acc in names(datasets)) {
  d<-datasets[[acc]];sc<-scores_for(d)
  write_tsv(score_effects(d,sc),paste0('04_results/sensitivity/',acc,'_common_members_module_effects.tsv'))
  write_tsv(data.frame(set=rownames(sc),sc,check.names=FALSE),paste0('04_results/sensitivity/',acc,'_common_members_module_scores.tsv'))
}
write_tsv(data.frame(set=names(allsets),n_common=lengths(allsets),genes=sapply(allsets,paste,collapse='|')),'04_results/sensitivity/common_module_members.tsv')
allsets<-original_allsets
saveRDS(list(data=datasets,scores=scorelist,sets=sets),'03_processed/bulk_analysis_objects.rds')
writeLines(capture.output(sessionInfo()),'env/R_sessionInfo.txt')
cat('COMPLETE\n')
