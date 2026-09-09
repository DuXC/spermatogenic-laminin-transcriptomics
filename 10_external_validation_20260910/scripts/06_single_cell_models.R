options(stringsAsFactors=FALSE,width=160)
.libPaths(c(normalizePath('../env/R_library'),.libPaths()))
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
set.seed(20260910)
stopifnot(file.exists('00_admin/external_annotation_freeze_receipt.json'))
dir.create('04_results/single_cell',recursive=TRUE,showWarnings=FALSE)
write_tsv<-function(x,name){p<-file.path('04_results/single_cell',name);if(grepl('\\.gz$',p)){con<-gzfile(p,'wt');on.exit(close(con));p<-con};write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
counts<-as.matrix(read.delim('03_processed/single_cell/all_pseudobulk_counts.tsv.gz',row.names=1,check.names=FALSE))
samples<-read.delim('03_processed/single_cell/all_pseudobulk_samples.tsv')
stopifnot(identical(colnames(counts),samples$sample_id),all(is.finite(counts)),all(counts>=0))
sets<-fromJSON('01_metadata/fixed_gene_sets.json')$primary
main<-c('Peritubular_compartment','Sertoli','Leydig','Endothelial','Macrophage','Perivascular')
coverage<-list();gene_tables<-list();camera_tables<-list();norm_qc<-list();descriptive_scores<-list();member_coverage<-list()
for(variant in c('primary','mt20'))for(comp in c(main,'Author_PMC_state','Author_fibrotic_PMC_state')){
  s<-samples[samples$variant==variant&samples$compartment==comp,]
  ok<-nrow(s)==6&&all(s$n_cells>=30)&&sum(s$group=='OA')==3&&sum(s$group=='Crypto')==3
  coverage[[length(coverage)+1]]<-data.frame(variant=variant,compartment=comp,n_OA_eligible=sum(s$group=='OA'&s$n_cells>=30),n_Crypto_eligible=sum(s$group=='Crypto'&s$n_cells>=30),minimum_cells=min(s$n_cells),eligible=ok)
  if(!ok||comp=='Author_fibrotic_PMC_state')next
  s$group<-factor(s$group,levels=c('OA','Crypto'))
  x<-counts[,s$sample_id,drop=FALSE];y<-DGEList(x);keep<-rowSums(cpm(y)>=1)>=3;y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM')
  for(age_adjusted in if(variant=='primary')c(FALSE,TRUE)else FALSE){
    design<-model.matrix(if(age_adjusted)~group+age else~group,s);stopifnot(qr(design)$rank==ncol(design))
    scenario<-paste(variant,comp,if(age_adjusted)'age_adjusted'else'unadjusted',sep='__')
    v<-voom(y,design,plot=FALSE);fit<-eBayes(lmFit(v,design),robust=TRUE)
    g<-topTable(fit,coef=2,n=Inf,sort.by='none',confint=TRUE);g$gene<-rownames(g);g$SE_moderated<-fit$stdev.unscaled[,2]*sqrt(fit$s2.post)
    g$scenario<-scenario;g$variant<-variant;g$compartment<-comp;g$age_adjusted<-age_adjusted;g$residual_df<-nrow(design)-ncol(design)
    gene_tables[[scenario]]<-g[,c('scenario','variant','compartment','age_adjusted','residual_df','gene','logFC','CI.L','CI.R','SE_moderated','AveExpr','t','P.Value','adj.P.Val')]
    idx<-ids2indices(sets,rownames(v));eligible_sets<-idx[lengths(idx)>=10]
    for(nm in names(sets))member_coverage[[length(member_coverage)+1]]<-data.frame(scenario=scenario,program=nm,n_total=length(sets[[nm]]),n_measured=length(idx[[nm]]),eligible=length(idx[[nm]])>=10,genes=paste(rownames(v)[idx[[nm]]],collapse='|'))
    for(mode in c('estimated','fixed_0.01')){
      ca<-camera(v,eligible_sets,design,contrast=2,inter.gene.cor=if(mode=='estimated')NA_real_ else .01,allow.neg.cor=FALSE,sort=FALSE)
      if(!'Correlation'%in%names(ca))ca$Correlation<-.01
      ca$program<-rownames(ca);ca$scenario<-scenario;ca$variant<-variant;ca$compartment<-comp;ca$age_adjusted<-age_adjusted;ca$mode<-mode
      ca$family<-if(variant=='primary'&&!age_adjusted&&comp%in%main)'primary_somatic'else scenario
      ca$mean_member_log2FC<-vapply(ca$program,function(nm)mean(g$logFC[idx[[nm]]]),numeric(1))
      ca$interpretation<-if(comp=='Peritubular_compartment')'author-defined mixed peritubular compartment'else'author PMC contractile-state sensitivity'
      camera_tables[[length(camera_tables)+1]]<-ca
    }
    if(!age_adjusted){
      write_tsv(data.frame(gene=rownames(v),v$E,check.names=FALSE),paste0(scenario,'_voom_logCPM.tsv.gz'))
      z<-t(scale(t(v$E)))
      for(nm in names(sets))if(length(idx[[nm]])>=3){
        val<-colMeans(z[idx[[nm]],,drop=FALSE]);descriptive_scores[[length(descriptive_scores)+1]]<-data.frame(scenario=scenario,variant=variant,compartment=comp,program=nm,donor=s$donor,group=s$group,n_cells=s$n_cells,score=val,n_members=length(idx[[nm]]))
      }
      norm_qc[[length(norm_qc)+1]]<-cbind(s,y$samples,genes_retained=nrow(y))
    }
    cat('COMPLETE',scenario,'genes',nrow(v),'residual_df',nrow(design)-ncol(design),'\n')
  }
}
write_tsv(do.call(rbind,coverage),'donor_compartment_eligibility.tsv')
g<-do.call(rbind,gene_tables);write_tsv(g,'all_gene_effects.tsv.gz')
ca<-do.call(rbind,camera_tables);ca$FDR_family<-NA_real_
for(family in unique(ca$family))for(mode in unique(ca$mode)){
  ii<-ca$family==family&ca$mode==mode;nt<-length(unique(ca$compartment[ii]))*6;ca$FDR_family[ii]<-p.adjust(ca$PValue[ii],'BH',n=nt)
}
write_tsv(ca,'program_camera_and_effects.tsv')
write_tsv(do.call(rbind,member_coverage),'program_member_coverage.tsv')
write_tsv(do.call(rbind,descriptive_scores),'donor_descriptive_program_scores.tsv')
write_tsv(do.call(rbind,norm_qc),'pseudobulk_library_QC.tsv')
lam<-sets[['REACTOME_LAMININ_INTERACTIONS']];out<-list()
for(scenario in names(gene_tables)){
  a<-gene_tables[[scenario]];x<-merge(data.frame(gene=lam),a,by='gene',all.x=TRUE,sort=FALSE);x$scenario<-scenario;x$measured<-!is.na(x$logFC);out[[scenario]]<-x
}
write_tsv(do.call(rbind,out),'all_30_laminin_gene_effects.tsv')
writeLines(capture.output(sessionInfo()),'00_admin/external_single_cell_R_sessionInfo.txt')
print(ca[ca$program=='REACTOME_LAMININ_INTERACTIONS'&ca$mode=='estimated',c('scenario','NGenes','Direction','PValue','FDR_family','mean_member_log2FC')])
cat('ALL_EXTERNAL_SINGLE_CELL_MODELS_COMPLETE\n')
