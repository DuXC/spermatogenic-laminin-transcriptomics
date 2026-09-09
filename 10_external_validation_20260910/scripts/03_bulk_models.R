options(stringsAsFactors=FALSE,width=160)
.libPaths(c(normalizePath('../env/R_library'),.libPaths()))
suppressPackageStartupMessages({library(limma);library(jsonlite)})
set.seed(20260910)
dir.create('04_results/bulk',recursive=TRUE,showWarnings=FALSE)
write_tsv<-function(x,name){p<-file.path('04_results/bulk',name);if(grepl('\\.gz$',p)){con<-gzfile(p,'wt');on.exit(close(con));p<-con};write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
read_expression<-function(acc,variant){as.matrix(read.delim(paste0('03_processed/',acc,'_',variant,'_gene_expression.tsv.gz'),row.names=1,check.names=FALSE))}
sets<-fromJSON('01_metadata/fixed_gene_sets.json');allsets<-c(sets$primary,sets$exploratory_markers)
primary_names<-names(sets$primary);main_effects<-list();main_camera<-list();statuses<-list()
score_matrix<-function(e){z<-t(scale(t(e)));ans<-lapply(allsets,function(g){g<-intersect(g,rownames(z));if(length(g)<3)rep(NA_real_,ncol(z))else colMeans(z[g,,drop=FALSE])});m<-do.call(rbind,ans);colnames(m)<-colnames(e);m}
score_fit<-function(sc,s,adjust=FALSE){
  design<-model.matrix(~group,s)
  if(adjust){if(anyNA(sc['Postmeiotic_markers',]))return(NULL);design<-cbind(design,postmeiotic=as.numeric(sc['Postmeiotic_markers',]))}
  if(qr(design)$rank<ncol(design))return(NULL)
  df<-nrow(design)-ncol(design);inv<-solve(crossprod(design));out<-list()
  for(nm in rownames(sc)){
    if(anyNA(sc[nm,])||(adjust&&nm=='Postmeiotic_markers'))next
    f<-lm.fit(design,sc[nm,]);se<-sqrt(sum(f$residuals^2)/df*inv[2,2]);est<-f$coefficients[2]
    out[[length(out)+1]]<-data.frame(cohort=s$cohort[1],program=nm,effect=est,SE=se,CI_low=est-qt(.975,df)*se,CI_high=est+qt(.975,df)*se,P=2*pt(-abs(est/se),df),df=df,postmeiotic_adjustment=adjust,n_NOA=sum(s$group=='NOA'),n_OA=sum(s$group=='OA'))
  }
  o<-do.call(rbind,out);o$FDR_six<-NA_real_;ii<-o$program%in%primary_names;o$FDR_six[ii]<-p.adjust(o$P[ii],'BH');o
}
genes<-function(e,s,variant){
  design<-model.matrix(~group,s);fit<-eBayes(lmFit(e,design),robust=TRUE,trend=TRUE)
  a<-topTable(fit,coef=2,n=Inf,sort.by='none',confint=TRUE);a$gene<-rownames(a)
  names(a)[names(a)=='logFC']<-'expression_difference';a$cohort<-s$cohort[1];a$variant<-variant;a$unit<-s$gene_effect_unit[1]
  a$SE_moderated<-fit$stdev.unscaled[,2]*sqrt(fit$s2.post)
  a$n_observed_OA<-rowSums(!is.na(e[,s$group=='OA',drop=FALSE]));a$n_observed_NOA<-rowSums(!is.na(e[,s$group=='NOA',drop=FALSE]))
  a[,c('cohort','variant','gene','unit','expression_difference','CI.L','CI.R','SE_moderated','AveExpr','t','P.Value','adj.P.Val','n_observed_OA','n_observed_NOA')]
}
for(acc in c('GSE9210','GSE108886')){
  s<-read.delim(paste0('01_metadata/',acc,'_analysis_samples.tsv'));s$group<-factor(s$group,levels=c('OA','NOA'))
  e<-read_expression(acc,'primary');stopifnot(identical(colnames(e),s$gsm),!anyNA(e))
  g<-genes(e,s,'primary');write_tsv(g,paste0(acc,'_all_genes.tsv.gz'))
  sc<-score_matrix(e);write_tsv(data.frame(program=rownames(sc),sc,check.names=FALSE),paste0(acc,'_scores.tsv'))
  ef<-score_fit(sc,s);main_effects[[acc]]<-ef;write_tsv(ef,paste0(acc,'_score_effects.tsv'))
  adjusted<-score_fit(sc,s,TRUE)
  if(is.null(adjusted))statuses[[length(statuses)+1]]<-data.frame(cohort=acc,analysis='postmeiotic_adjustment',status='NOT_EVALUABLE',reason='Fewer than three measured postmeiotic markers or rank-deficient design') else write_tsv(adjusted,paste0(acc,'_postmeiotic_adjusted_scores.tsv'))
  index<-ids2indices(sets$primary,rownames(e));design<-model.matrix(~group,s)
  ca_all<-list()
  for(mode in c('estimated','fixed_0.01')){
    eligible<-index[lengths(index)>=10]
    ca<-camera(e,eligible,design,contrast=2,inter.gene.cor=if(mode=='estimated')NA_real_ else .01,allow.neg.cor=FALSE,trend.var=TRUE,sort=FALSE)
    if(!'Correlation'%in%names(ca))ca$Correlation<-.01
    ca$program<-rownames(ca);ca$cohort<-acc;ca$mode<-mode;ca$status<-'TESTED';ca$FDR_six<-p.adjust(ca$PValue,'BH',n=6)
    ca_all[[mode]]<-ca
    for(nm in names(index)[lengths(index)<10])statuses[[length(statuses)+1]]<-data.frame(cohort=acc,analysis=paste('camera',mode,nm,sep=':'),status='NOT_EVALUABLE',reason=paste('Measured members',length(index[[nm]]),'below original minimum 10'))
  }
  ca<-do.call(rbind,ca_all);main_camera[[acc]]<-ca;write_tsv(ca,paste0(acc,'_camera.tsv'))
  influence<-lapply(seq_len(nrow(s)),function(i){x<-score_fit(sc[,-i,drop=FALSE],s[-i,]);x$omitted_gsm<-s$gsm[i];x})
  write_tsv(do.call(rbind,influence),paste0(acc,'_leave_one_out.tsv'))
  variant<-if(acc=='GSE9210')'available_case'else'all_mapped';es<-read_expression(acc,variant)
  write_tsv(genes(es,s,variant),paste0(acc,'_',variant,'_genes.tsv.gz'))
  if(!anyNA(es))write_tsv(score_fit(score_matrix(es),s),paste0(acc,'_',variant,'_score_effects.tsv'))
  pc<-prcomp(t(e),center=TRUE,scale.=FALSE);write_tsv(cbind(s,pc$x[,1:3]),paste0(acc,'_PCA.tsv'))
  write_tsv(data.frame(gsm=s$gsm,cor(e),check.names=FALSE),paste0(acc,'_sample_correlations.tsv'))
  cat('COMPLETE',acc,'genes',nrow(e),'\n');print(ef[ef$program=='REACTOME_LAMININ_INTERACTIONS',c('program','effect','CI_low','CI_high','P','FDR_six')])
}
ef<-do.call(rbind,main_effects);ef$FDR_external_12<-NA_real_;ii<-ef$program%in%primary_names;ef$FDR_external_12[ii]<-p.adjust(ef$P[ii],'BH',n=12)
write_tsv(ef,'external_score_effects.tsv')
ca<-do.call(rbind,main_camera);ca$FDR_external_12<-NA_real_
for(mode in unique(ca$mode)){ii<-ca$mode==mode;ca$FDR_external_12[ii]<-p.adjust(ca$PValue[ii],'BH',n=12)}
write_tsv(ca,'external_camera.tsv');write_tsv(do.call(rbind,statuses),'analysis_availability.tsv')
laminin<-sets$primary[['REACTOME_LAMININ_INTERACTIONS']];rows<-list()
for(acc in c('GSE9210','GSE108886'))for(variant in c('primary',if(acc=='GSE9210')'available_case'else'all_mapped')){
  path<-paste0('04_results/bulk/',acc,if(variant=='primary')'_all_genes.tsv.gz'else paste0('_',variant,'_genes.tsv.gz'))
  a<-read.delim(path);x<-merge(data.frame(gene=laminin),a,by='gene',all.x=TRUE,sort=FALSE);x$cohort<-acc;x$variant<-variant;x$measured<-!is.na(x$expression_difference);rows[[length(rows)+1]]<-x
}
write_tsv(do.call(rbind,rows),'all_30_laminin_gene_evidence.tsv')
writeLines(capture.output(sessionInfo()),'00_admin/external_bulk_R_sessionInfo.txt')
cat('ALL_EXTERNAL_BULK_MODELS_COMPLETE\n')
