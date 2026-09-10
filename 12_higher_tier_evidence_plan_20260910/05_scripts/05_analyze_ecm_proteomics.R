args <- commandArgs(trailingOnly=TRUE)
N <- normalizePath(args[1])
P <- dirname(N)
.libPaths(c(file.path(P,"env/R_library"),.libPaths()))
suppressPackageStartupMessages(library(limma))
O <- file.path(N,"03_analysis/PXD011817")
readx <- function(f) as.matrix(read.delim(file.path(O,f),row.names=1,check.names=FALSE))
x <- readx("gene_log2_LFQ_min3_each_group.tsv")
xc <- readx("gene_log2_LFQ_complete_10.tsv")
samples <- read.delim(file.path(O,"participants.tsv"))
stopifnot(identical(colnames(x),samples$participant),all(is.finite(xc)),ncol(xc)==10)
g <- as.numeric(samples$group=="SRneg")
design <- cbind(Intercept=1,SRneg_minus_SRpos=g)
coverage <- read.delim(file.path(O,"fixed_program_coverage.tsv"),check.names=FALSE)
sets <- setNames(strsplit(coverage$genes,";",fixed=TRUE),coverage$program)
idx <- lapply(sets,function(v)which(rownames(xc)%in%v))
evaluable <- lengths(idx)>=10
camera_once <- function(xx,dd,cor) {
  cc <- camera(xx,index=idx[evaluable],design=dd,contrast=2,inter.gene.cor=cor)
  out <- data.frame(program=names(idx),members=lengths(idx),status=ifelse(evaluable,"EVALUABLE","NOT_EVALUABLE_COVERAGE"),Direction=NA_character_,PValue=NA_real_,FDR_six=NA_real_)
  out$Direction[evaluable] <- cc[out$program[evaluable],"Direction"]
  out$PValue[evaluable] <- cc[out$program[evaluable],"PValue"]
  family_p <- out$PValue;family_p[!evaluable] <- 1
  out$FDR_six <- p.adjust(family_p,"BH");out$FDR_six[!evaluable] <- NA_real_
  out
}
primary <- camera_once(xc,design,NA_real_)
fixed <- camera_once(xc,design,.01)
write.table(primary,file.path(O,"camera_estimated_correlation.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
write.table(fixed,file.path(O,"camera_fixed_correlation_001.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
fit <- eBayes(lmFit(x,design),robust=TRUE,trend=TRUE)
tt <- topTable(fit,coef=2,n=Inf,sort.by="none",confint=.95)
tt$gene <- rownames(tt)
tt$SRpos_n <- rowSums(!is.na(x[,g==0,drop=FALSE]))
tt$SRneg_n <- rowSums(!is.na(x[,g==1,drop=FALSE]))
write.table(tt,file.path(O,"all_gene_results.tsv"),sep="\t",row.names=FALSE,quote=FALSE,na="NA")
z <- t(scale(t(xc)))
stopifnot(all(is.finite(z)))
scoremat <- sapply(idx,function(v)colMeans(z[v,,drop=FALSE]))
scorefit <- function(y,gg) {
  fit <- lm(y~gg);co <- coef(summary(fit));ci <- confint(fit)
  c(effect=co[2,1],CI_low=ci[2,1],CI_high=ci[2,2],PValue=co[2,4])
}
scores <- as.data.frame(t(sapply(seq_along(idx),function(i) scorefit(scoremat[,i],g))))
scores$program <- names(idx);scores$members <- lengths(idx);scores$FDR_six <- p.adjust(scores$PValue,"BH")
write.table(scores,file.path(O,"program_score_group_effects.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(cbind(samples,scoremat),file.path(O,"participant_program_scores.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
lam <- "REACTOME_LAMININ_INTERACTIONS"
loo <- lapply(seq_len(ncol(xc)),function(i) {
  score <- scorefit(scoremat[-i,lam],g[-i])
  cc <- camera_once(xc[,-i,drop=FALSE],design[-i,,drop=FALSE],NA_real_)
  cc <- cc[cc$program==lam,]
  data.frame(omitted=samples$participant[i],group=samples$group[i],t(score),camera_direction=cc$Direction,camera_p=cc$PValue,camera_FDR_six=cc$FDR_six)
})
write.table(do.call(rbind,loo),file.path(O,"laminin_leave_one_participant_out.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(O,"R_sessionInfo.txt"))
print(primary)
print(scores[,c("program","members","effect","CI_low","CI_high","PValue","FDR_six")])
print(tt[tt$gene%in%sets[[lam]],c("gene","logFC","CI.L","CI.R","P.Value","adj.P.Val")])
