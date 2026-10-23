####################
# Analysis registry:
# Status: active exploratory terminal validation using existing bulk DNA profiles.
# Script: analysis/cnv/Auto_PDO_wes_cna_mp_validation.R
# Methodology: analysis/methodology/cnv/Auto_PDO_wes_cna_mp_validation_methodology.md
# Map: analysis/ANALYSIS_MAP.md
# Description: Compare six matched untreated PDO bulk-WES gene/arm CNA profiles
#   with sample-mean UCell activity of the exact 17 scATLAS MPs. Add PDO panels
#   to the fixed existing TCGA/scATLAS gene and arm comparisons.
# Inputs: PDOs_outs/PDOs_list_PDOs.rds; external sarek_mutect/variant_calling/
#   cnvkit/PDO_<id>_vs_NT_<id>/PDO_<id>.somatic.call.cns for IDs 1070,1072,
#   1090,1121,1141,1181; ITH_all/all_samples/hg38_gencode_v27.txt;
#   scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/
#   intermediate/merged_refined_mp_genes.rds and tables/centred_refined_mp_state_grouping.csv;
#   scRef_Pipeline/ref_outs/TCGA/cna_gene_event_validation/intermediate/Auto_targeted_cna_results.rds;
#   scRef_Pipeline/ref_outs/TCGA/cna_recurrent_event_validation/intermediate/Auto_tcga_cna_plot_source_data.rds;
#   scRef_Pipeline/ref_outs/TCGA/cna_gene_event_validation/intermediate/source_data/software/pymupdf_1.26.7/.
# Outputs: live PDOs_outs/Auto_PDO_wes_cna_mp_validation/{intermediate,tables,
#   figures,logs,reports}/Auto_* (scores, measurements, tests, figure sources,
#   figure cache, five PDFs, PDF audit and previews).
# Downstream: terminal exploratory validation; does not alter states or CNA calls.
# Cache: PDO_FORCE_REBUILD=1 rebuilds scores; PDO_REPLOT_ONLY=1 uses saved
#   plot objects and reruns PDF validation. No source CNA correction is applied.
# Run: /opt/pbs/bin/qsub job_scripts/cnv/Auto_PDO_wes_cna_mp_validation.sh
# Environment/resources: gnmf; 2 CPUs, 64 GB, 02:00:00; project-root PBS logs.
####################
library(Seurat)
library(UCell)
library(Matrix)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
source("analysis/shared/Auto_pdo_analysis_config.R")
source("analysis/shared/Auto_pdo_analysis_helpers.R")
if (!nzchar(Sys.getenv("PBS_JOBID"))) stop("Run on a PBS compute node.")
set.seed(42)
script_path <- "analysis/cnv/Auto_PDO_wes_cna_mp_validation.R"
out_dir <- file.path(PDO_OUTPUT_DIR, "Auto_PDO_wes_cna_mp_validation")
tiers <- pdo_ensure_output_tiers(out_dir)
figure_cache <- file.path(tiers[["intermediate"]], "Auto_figure_bundle.rds")
audit_pdf <- function() {
  python_bin <- file.path(Sys.getenv("CONDA_PREFIX"),"bin","python3")
  if (!nzchar(Sys.getenv("CONDA_PREFIX")) || !file.exists(python_bin)) stop("Active analysis environment Python is required.")
  status <- system2(python_bin, c("analysis/cnv/Auto_validate_PDO_wes_cna_figures.py", shQuote(out_dir)))
  if (status != 0L) stop("PDF audit failed: ", status)
}
save_figures <- function(bundle) {
  for (item in bundle) {
    path <- file.path(tiers[["figures"]], item$name)
    grDevices::cairo_pdf(path, width=item$width, height=item$height, onefile=TRUE)
    for (page in item$pages) {
      grid::grid.newpage()
      if (is.function(page)) page <- page()
      if (inherits(page,"ggplot")) print(page,newpage=FALSE) else grid::grid.draw(page)
    }
    grDevices::dev.off()
  }
}
if (pdo_get_env_flag("PDO_REPLOT_ONLY")) {
  if (pdo_get_env_flag("PDO_FORCE_REBUILD")) stop("Rebuild and replot flags conflict.")
  pdo_require_files(figure_cache)
  save_figures(readRDS(figure_cache)); audit_pdf()
  message("Replot and PDF audit passed: ", out_dir)
  quit(save="no")
}
ids <- c("1070", "1072", "1090", "1121", "1141", "1181")
samples <- paste0("SUR", ids, "_Untreated_PDO")
sarek <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect"
sc_root <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs"
cn_paths <- file.path(sarek,"variant_calling/cnvkit",paste0("PDO_",ids,"_vs_NT_",ids),paste0("PDO_",ids,".somatic.call.cns"))
mp_path <- file.path(sc_root,"Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds")
group_path <- file.path(sc_root,"Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv")
reference_path <- file.path(sc_root,"TCGA/cna_gene_event_validation/intermediate/Auto_targeted_cna_results.rds")
parent_path <- file.path(sc_root,"TCGA/cna_recurrent_event_validation/intermediate/Auto_tcga_cna_plot_source_data.rds")
gene_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
rna_path <- file.path(PDO_OUTPUT_DIR,"PDOs_list_PDOs.rds")
inputs <- c(cn_paths,mp_path,group_path,reference_path,parent_path,gene_path,rna_path)
pdo_require_files(inputs)
reference <- readRDS(reference_path)
parent <- readRDS(parent_path)
mp_sets <- readRDS(mp_path)
mp_order <- names(reference$mp_sets)
grouping <- read.csv(group_path)
stopifnot(length(mp_order)==17L, setequal(mp_order,grouping$mp),
          !any(c("MP2x","MP11c","MP18a") %in% mp_order),
          setequal(names(mp_sets),mp_order),identical(mp_sets[mp_order], reference$mp_sets))
mp_sets <- mp_sets[mp_order]
mp_labels <- setNames(grouping$plot_label,grouping$mp)
genes <- reference$selected_genes
stopifnot(nrow(genes)==20L,all(table(genes$direction)==10L),nrow(reference$top_gene)==10L,nrow(reference$top_arm)==10L)
manifest <- data.frame(donor=paste0("SUR",ids),sample=samples,dna_path=cn_paths,
  pairing="User confirmed same sample by concordant WES/RNA names; untreated RNA used once per donor",
  dna_unit="bulk PDO culture",rna_unit="sample mean of single-cell UCell scores",
  limitation="Existing CNVkit WGS-mode output from WES; not reprocessed; exploratory",stringsAsFactors=FALSE)
write.csv(manifest,file.path(tiers[["tables"]],"Auto_sample_manifest.csv"),row.names=FALSE)
write.csv(genes,file.path(tiers[["tables"]],"Auto_selected_oac_genes.csv"),row.names=FALSE)
write.csv(grouping,file.path(tiers[["tables"]],"Auto_scatlas_mp_definitions.csv"),row.names=FALSE)

# Score all six RNA samples with one prespecified method, retaining exact cell IDs.
score_path <- file.path(tiers[["intermediate"]],"Auto_scatlas_mp_cell_scores.rds")
fingerprint <- data.frame(path=c(rna_path,mp_path,group_path),size=file.info(c(rna_path,mp_path,group_path))$size,
  mtime=as.numeric(file.info(c(rna_path,mp_path,group_path))$mtime))
if (file.exists(score_path) && !pdo_get_env_flag("PDO_FORCE_REBUILD")) {
  score_cache <- readRDS(score_path)
  if (!identical(score_cache$fingerprint,fingerprint) || !identical(score_cache$samples,samples)) stop("Score cache inputs changed; explicitly rebuild.")
} else {
  pdos <- readRDS(rna_path)
  if (!all(samples %in% names(pdos))) stop("Missing RNA samples.")
  pdos <- pdos[samples]
  score_rows <- list(); coverage_rows <- list()
  for (sample_id in samples) {
    message("UCell: ",sample_id)
    obj <- pdos[[sample_id]]
    if (!all(as.character(obj@meta.data$orig.ident)==sample_id)) stop("orig.ident mismatch: ",sample_id)
    counts <- Seurat::GetAssayData(obj,assay="RNA",layer="counts")
    stopifnot(ncol(counts)>0L,!anyDuplicated(colnames(counts)),all(Matrix::colSums(counts)>0))
    features <- lapply(mp_sets,intersect,rownames(counts))
    if (any(lengths(features)<5L)) stop("Insufficient MP gene coverage: ",sample_id)
    coverage_rows[[sample_id]] <- data.frame(sample=sample_id,mp=mp_order,n_signature=lengths(mp_sets),n_measured=lengths(features))
    scores <- as.data.frame(UCell::ScoreSignatures_UCell(matrix=counts,features=features,maxRank=1500,
      chunk.size=1000,ncores=2,force.gc=TRUE))
    names(scores) <- sub("_UCell$","",names(scores))
    if (!setequal(names(scores),mp_order) || !setequal(rownames(scores),colnames(counts))) stop("UCell identity mismatch.")
    score_rows[[sample_id]] <- data.frame(cell=rownames(scores),sample=sample_id,scores[,mp_order,drop=FALSE],check.names=FALSE)
    pdos[[sample_id]] <- NULL
    rm(obj,counts,scores); gc(verbose=FALSE)
  }
  score_cache <- list(fingerprint=fingerprint,samples=samples,cells=bind_rows(score_rows),coverage=bind_rows(coverage_rows),
    method="UCell maxRank=1500, full measured RNA feature universe per sample",session=sessionInfo())
  saveRDS(score_cache,score_path)
  rm(pdos); gc(verbose=FALSE)
}
cells <- score_cache$cells
stopifnot(all(is.finite(as.matrix(cells[,mp_order]))),!anyDuplicated(cells[c("sample","cell")]))
scores <- cells |> group_by(sample) |> summarise(n_cells=n(),across(all_of(mp_order),mean),.groups="drop")
scores <- scores[match(samples,scores$sample),]
stopifnot(identical(scores$sample,samples),nrow(scores)==6L)
write.csv(scores,file.path(tiers[["tables"]],"Auto_sample_mp_activity.csv"),row.names=FALSE)
write.csv(score_cache$coverage,file.path(tiers[["tables"]],"Auto_mp_gene_coverage.csv"),row.names=FALSE)

# Read original DNA segments. No RNA-dependent centering or absolute-CN shift.
segments <- bind_rows(lapply(seq_along(samples),function(i) {
  x <- fread(cn_paths[i],data.table=FALSE)
  stopifnot(all(c("chromosome","start","end","log2","cn","probes") %in% names(x)),
    all(is.finite(x$log2)),all(x$end>x$start),all(x$cn>=0))
  x$sample <- samples[i]; x
}))
go <- read.table(gene_path,header=FALSE,col.names=c("gene","chr","start","end"),stringsAsFactors=FALSE)
go <- go[!duplicated(go$gene),]
if (!all(genes$gene %in% go$gene)) stop("Selected gene coordinate missing.")
gene_values <- bind_rows(lapply(samples,function(sample_id) bind_rows(lapply(genes$gene,function(gene_id) {
  g <- go[go$gene==gene_id,]; s <- segments[segments$sample==sample_id & segments$chromosome==g$chr,]
  # Gene order is 1-based inclusive; CNVkit segments are 0-based half-open.
  overlap <- pmax(0,pmin(s$end,g$end)-pmax(s$start,g$start-1))
  hit <- overlap>0
  covered <- sum(overlap[hit])/(g$end-g$start+1)
  if (covered>1+1e-8) stop("Overlapping DNA segments at ",sample_id," ",gene_id)
  data.frame(sample=sample_id,gene=gene_id,chromosome=g$chr,gene_start=g$start,gene_end=g$end,
    covered_fraction=covered,n_segments=sum(hit),n_segment_probes=sum(s$probes[hit]),
    value=if(any(hit)) weighted.mean(s$log2[hit],overlap[hit]) else NA_real_,
    cn=if(sum(hit)==1L && covered>=0.999) s$cn[hit] else NA_real_)
}))))
if (any(gene_values$covered_fraction<0.9)) stop("A selected gene has less than 90% segment coverage; do not impute neutral.")
centromeres <- c(chr1=121700000,chr2=91800000,chr3=87900000,chr4=50600000,chr5=48400000,chr6=61000000,
 chr7=59900000,chr8=45600000,chr9=49000000,chr10=40200000,chr11=53400000,chr12=35500000,
 chr13=17700000,chr14=17200000,chr15=19000000,chr16=36800000,chr17=25100000,chr18=18500000,
 chr19=26200000,chr20=28100000,chr21=12000000,chr22=15000000,chrX=61000000)
sizes <- c(chr1=248956422,chr2=242193529,chr3=198295559,chr4=190214555,chr5=181538259,chr6=170805979,
 chr7=159345973,chr8=145138636,chr9=138394717,chr10=133797422,chr11=135086622,chr12=133275309,
 chr13=114364328,chr14=107043718,chr15=101991189,chr16=90338345,chr17=83257441,chr18=80373285,
 chr19=58617616,chr20=64444167,chr21=46709983,chr22=50818468,chrX=156040895)
arms <- as.vector(rbind(paste0(names(sizes),"p"),paste0(names(sizes),"q")))
arm_values <- bind_rows(lapply(samples,function(sample_id) bind_rows(lapply(arms,function(arm_id) {
  chr <- sub("[pq]$","",arm_id); left <- if(endsWith(arm_id,"p")) 0 else centromeres[chr]
  right <- if(endsWith(arm_id,"p")) centromeres[chr] else sizes[chr]
  s <- segments[segments$sample==sample_id & segments$chromosome==chr,]
  overlap <- pmax(0,pmin(s$end,right)-pmax(s$start,left)); hit <- overlap>0
  data.frame(sample=sample_id,arm=arm_id,value=if(any(hit)) weighted.mean(s$log2[hit],overlap[hit]) else NA_real_,
    covered_fraction=sum(overlap)/(right-left),n_segments=sum(hit))
}))))
if (any(!is.finite(arm_values$value))) stop("Unmeasured chromosome arm; no neutral imputation.")
write.csv(gene_values,file.path(tiers[["tables"]],"Auto_gene_cna_measurements.csv"),row.names=FALSE)
write.csv(arm_values,file.path(tiers[["tables"]],"Auto_arm_cna_measurements.csv"),row.names=FALSE)
write.csv(segments,file.path(tiers[["tables"]],"Auto_original_dna_segments.csv"),row.names=FALSE)

# Exhaustive two-sided permutation P values: six independent PDOs, including ties.
permutations <- function(v) {
  if(length(v)==1L) return(matrix(v,nrow=1L))
  do.call(rbind,lapply(seq_along(v),function(i)cbind(v[i],permutations(v[-i]))))
}
perms <- permutations(seq_len(6L))
stopifnot(nrow(perms)==720L)
exact_test <- function(x,y) {
  stopifnot(length(x)==6L,length(y)==6L,all(is.finite(x)),all(is.finite(y)))
  if(length(unique(x))<2L || length(unique(y))<2L) return(data.frame(n=6L,rho=NA_real_,p_value=NA_real_,n_permutations=720L))
  rx <- rank(x); ry <- rank(y); rho <- cor(rx,ry)
  null <- apply(perms,1,function(p)cor(rx,ry[p]))
  data.frame(n=6L,rho=rho,p_value=mean(abs(null)>=abs(rho)-1e-12),n_permutations=720L)
}
make_tests <- function(values,key,ref_tests) {
  tests <- bind_rows(lapply(unique(values[[key]]),function(feature) {
    x <- values[values[[key]]==feature,]; x <- x[match(samples,x$sample),]
    bind_rows(lapply(mp_order,function(mp)cbind(data.frame(feature=feature,mp=mp),exact_test(x$value,scores[[mp]]))))
  }))
  tests$fdr <- p.adjust(tests$p_value,"BH",n=length(unique(values[[key]]))*17L)
  names(tests)[names(tests)=="feature"] <- key
  tests$dataset <- "PDO bulk WES CNA"
  tests$predictor <- "Original segment log2 ratio; unshifted"
  ref_tests <- as.data.frame(ref_tests)
  ref_tests$n_permutations <- NA_integer_; ref_tests$predictor <- "Pinned reference measurement"
  bind_rows(ref_tests,tests)
}
gene_tests <- make_tests(gene_values,"gene",reference$gene_tests)
arm_tests <- make_tests(arm_values,"arm",reference$arm_tests)
write.csv(gene_tests,file.path(tiers[["tables"]],"Auto_gene_mp_tests.csv"),row.names=FALSE)
write.csv(arm_tests,file.path(tiers[["tables"]],"Auto_arm_mp_tests.csv"),row.names=FALSE)

# Three-cohort summary rectangles retain reference geometry and TCGA anchor.
rectangle_sources <- list()
tile_page <- function(tests,key,features,title,view) {
  df <- tests[tests[[key]] %in% features,]
  anchor <- df[grepl("^TCGA",df$dataset),c(key,"mp","rho","fdr")]
  names(anchor)[3:4] <- c("tcga_rho","tcga_fdr")
  df <- left_join(df,anchor,by=c(key,"mp"))
  df$same_direction <- is.finite(df$rho)&is.finite(df$tcga_rho)&df$rho*df$tcga_rho>0
  df$tcga_significant <- is.finite(df$tcga_fdr)&df$tcga_fdr<0.05
  df$view <- view
  df$display_rho <- ifelse(view=="All trends"|df$tcga_significant,df$rho,NA_real_)
  df$symbol <- ifelse(df$tcga_significant&df$same_direction,"*","")
  if(key=="gene") {
    sc_label <- "scATLAS gene RNA expression"
    axis_labels <- setNames(genes$event_label,genes$gene)
  } else { sc_label <- "scATLAS manual InferCNA"; axis_labels <- setNames(features,features) }
  tc_label <- unique(df$dataset[grepl("^TCGA",df$dataset)])
  df$dataset <- factor(df$dataset,levels=c(sc_label,tc_label,"PDO bulk WES CNA"))
  df$feature <- factor(df[[key]],levels=features)
  df$mp_label <- factor(mp_labels[df$mp],levels=rev(mp_labels[mp_order]))
  rectangle_sources[[length(rectangle_sources)+1L]] <<- df
  p <- ggplot(df,aes(feature,mp_label))+geom_tile(aes(fill=display_rho),colour="grey85",linewidth=0.4)+
    geom_text(aes(label=symbol),size=6)+facet_grid(.~dataset)+
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-0.25,0.25),
      oob=scales::squish,na.value="white",name="Spearman rho")+
    scale_x_discrete(labels=axis_labels,drop=FALSE)+
    labs(title=paste(title,view,sep=": "),x=NULL,y=NULL,
      subtitle="PDO: six matched untreated samples; original WES segments, exploratory",
      caption="* TCGA BH FDR < 0.05 and matching nonzero direction; does not imply PDO significance. White: masked/undefined. Colour saturates at +/-0.25.")+
    theme_classic(base_size=16)+theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=if(key=="arm") 9 else 13),
      axis.text.y=element_text(colour="black",size=12),strip.text=element_text(face="bold",size=16),
      plot.title=element_text(face="bold",size=22),plot.subtitle=element_text(size=15),plot.caption=element_text(size=12),
      panel.border=element_rect(colour="black",fill=NA),axis.line=element_blank(),plot.margin=margin(16,24,16,16))
  p
}
gene_pages <- list()
for(direction in c("gain","loss")) for(view in c("All trends","TCGA significant only")) {
  features <- genes$gene[as.character(genes$direction)==direction]
  gene_pages[[length(gene_pages)+1L]] <- tile_page(gene_tests,"gene",features,paste("Recurrent OAC",direction,"genes versus MP activity"),view)
}
arm_pages <- lapply(c("All trends","TCGA significant only"),function(view)tile_page(arm_tests,"arm",arms,"Chromosome-arm CNA versus MP activity",view))
write.csv(bind_rows(rectangle_sources),file.path(tiers[["tables"]],"Auto_rectangle_source.csv"),row.names=FALSE)

# Fixed top ten pairs from the existing TCGA/scATLAS run, never selected on PDO results.
box_sources <- list(); comparison_rows <- list(); pair_summaries <- list(); scatter_pages <- list()
group_colours <- c("Deep loss"="#053061",Loss="#2166AC",Neutral="grey72",Gain="#B2182B",Amplified="#67001F",Low="#2166AC",Medium="grey72",High="#B2182B")
stat_text <- function(row,pdo=FALSE) paste0("rho = ",sprintf("%.2f",row$rho),"; ",if(pdo) "exact P" else "P"," = ",
  format.pval(row$p_value,digits=2),"; BH FDR = ",format.pval(row$fdr,digits=2),"; n = ",row$n)
box_panel <- function(df,title,statistic,pdo=FALSE,x_label=NULL) {
  counts <- table(df$group)
  p <- ggplot(df,aes(group,y,fill=group))+geom_boxplot(outlier.shape=NA,width=0.62,linewidth=0.8,alpha=0.88)+
    geom_point(position=position_jitter(width=0.10,height=0,seed=42),size=if(pdo) 3.5 else 2.1,alpha=if(pdo) 0.9 else 0.4)+
    scale_fill_manual(values=group_colours,drop=FALSE)+scale_x_discrete(drop=FALSE,
      labels=setNames(paste0(names(counts),"\n(n=",counts,")"),names(counts)))+
    labs(title=title,subtitle=statistic,x=x_label,y=if(grepl("TCGA",title)) "Bulk GSVA MP activity" else
      if(grepl("scATLAS",title)) "Subclone MP activity (reference scale)" else "Mean UCell MP activity")+
    theme_classic(base_size=19)+theme(legend.position="none",plot.title=element_text(face="bold",size=21),
      plot.subtitle=element_text(size=13),axis.text=element_text(size=16,colour="black"),plot.margin=margin(14,18,14,18))
  if(pdo) p <- p+geom_text(aes(label=sub("_Untreated_PDO$","",sample)),
    position=position_jitter(width=0.10,height=0,seed=42),vjust=-0.8,size=3.2,check_overlap=FALSE)+
    scale_y_continuous(expand=expansion(mult=c(0.08,0.22)))
  p
}
make_triptych <- function(left,middle,right,title) {
  force(left); force(middle); force(right); force(title)
  function() arrangeGrob(left,middle,right,ncol=3,
    top=textGrob(title,gp=gpar(fontsize=20,fontface="bold")),
    bottom=textGrob("Fixed TCGA/scATLAS-selected pair; six independent PDOs. Original WES segments are exploratory. PDO P tests continuous log2 CNA.",gp=gpar(fontsize=13)))
}
pair_pages <- function(kind) {
  selected <- if(kind=="gene") reference$top_gene else reference$top_arm
  tests <- if(kind=="gene") gene_tests else arm_tests
  values <- if(kind=="gene") gene_values else arm_values
  result <- list()
  for(i in seq_len(nrow(selected))) {
    row <- selected[i,]; feature <- row[[kind]]; mp <- row$mp
    if(kind=="gene") {
      sc <- reference$sc_gene[reference$sc_gene$gene==feature,]
      sc <- data.frame(sample=sc$sample,x=sc$expression,y=sc[[paste0("mp__",mp)]])
      cuts <- quantile(sc$x,c(1/3,2/3),type=7,na.rm=TRUE)
      sc$group <- factor(ifelse(sc$x<=cuts[1],"Low",ifelse(sc$x<=cuts[2],"Medium","High")),levels=c("Low","Medium","High"))
      tc <- reference$tcga_gene[reference$tcga_gene$gene==feature,]
      tc <- data.frame(sample=tc$sample_barcode,x=tc$cna,y=tc[[mp]])
      tc$group <- factor(tc$x,levels=-2:2,labels=c("Deep loss","Loss","Neutral","Gain","Amplified"))
    } else {
      sc <- reference$sc_arm
      sc <- data.frame(sample=sc$sample,x=sc[[feature]],y=sc[[paste0("mp__",mp)]])
      sc$group <- factor(ifelse(sc$x>=0.05,"Gain",ifelse(sc$x<=-0.05,"Loss","Neutral")),levels=c("Loss","Neutral","Gain"))
      tc <- reference$tcga_arm[reference$tcga_arm$arm==feature,]
      tc <- data.frame(sample=tc$sample_barcode,x=tc$arm_mean,y=tc[[mp]])
      tc$group <- factor(ifelse(tc$x>=parent$selected_threshold,"Gain",ifelse(tc$x<=-parent$selected_threshold,"Loss","Neutral")),levels=c("Loss","Neutral","Gain"))
    }
    pdo <- values[values[[kind]]==feature,]; pdo <- pdo[match(samples,pdo$sample),]
    pdo$x <- pdo$value; pdo$y <- scores[[mp]]
    if(kind=="gene") pdo$group <- factor(ifelse(is.na(pdo$cn),NA,ifelse(pdo$cn>=4,"Amplified",
      c("Deep loss","Loss","Neutral","Gain")[pdo$cn+1])),levels=c("Deep loss","Loss","Neutral","Gain","Amplified"))
    else pdo$group <- factor(ifelse(pdo$x>=parent$selected_threshold,"Gain",ifelse(pdo$x<=-parent$selected_threshold,"Loss","Neutral")),levels=c("Loss","Neutral","Gain"))
    if(anyNA(pdo$group)) stop("A plotted gene spans incompatible CNA categories: ",feature)
    df_tests <- tests[tests[[kind]]==feature & tests$mp==mp,]
    sc_stat <- df_tests[grepl("^scATLAS",df_tests$dataset),]
    tc_stat <- df_tests[grepl("^TCGA",df_tests$dataset),]
    pdo_stat <- df_tests[df_tests$dataset=="PDO bulk WES CNA",]
    stopifnot(nrow(sc_stat)==1L,nrow(tc_stat)==1L,nrow(pdo_stat)==1L,nrow(pdo)==6L)
    pair_summaries[[paste(kind,i)]] <<- data.frame(kind=kind,feature=feature,mp=mp,reference_rank=i,
      pdo_rho=pdo_stat$rho,pdo_p=pdo_stat$p_value,pdo_fdr=pdo_stat$fdr,tcga_rho=tc_stat$rho,
      same_direction_tcga=is.finite(pdo_stat$rho)&&is.finite(tc_stat$rho)&&pdo_stat$rho*tc_stat$rho>0,
      gene_in_mp=kind=="gene" && feature %in% mp_sets[[mp]])
    cohort <- list(scATLAS=sc,TCGA=tc,PDO=pdo)
    for(label in names(cohort)) {
      d <- cohort[[label]]
      box_sources[[paste(kind,i,label)]] <<- data.frame(kind=kind,feature=feature,mp=mp,dataset=label,
        sample=d$sample,x=d$x,y=d$y,group=as.character(d$group))
      present <- names(table(d$group))[table(d$group)>=2L]
      if(length(present)>1L) for(pair in combn(present,2,simplify=FALSE)) {
        a <- d$y[d$group==pair[1]]; b <- d$y[d$group==pair[2]]
        p <- if(length(unique(c(a,b)))>1L) suppressWarnings(wilcox.test(a,b,exact=FALSE)$p.value) else NA_real_
        comparison_rows[[length(comparison_rows)+1L]] <<- data.frame(kind=kind,feature=feature,mp=mp,dataset=label,
          group1=pair[1],group2=pair[2],n1=length(a),n2=length(b),p_value=p)
      }
    }
    title <- paste0(i,". ",feature," versus ",mp_labels[mp],
      if(kind=="gene" && feature %in% mp_sets[[mp]]) " [gene belongs to MP]" else "")
    left <- box_panel(sc,if(kind=="gene") "scATLAS: gene RNA" else "scATLAS: inferred arm CNA",stat_text(sc_stat),
      x_label=if(kind=="gene") "RNA expression tertiles" else "InferCNA arm category")
    middle <- box_panel(tc,if(kind=="gene") "TCGA: gene CNA" else "TCGA: arm CNA",stat_text(tc_stat),
      x_label=if(kind=="gene") "GISTIC category" else paste0("Arm log2 threshold +/-",parent$selected_threshold))
    right <- box_panel(pdo,"PDO: bulk WES CNA",stat_text(pdo_stat,TRUE),TRUE,
      x_label=if(kind=="gene") "CNVkit estimated CN: 0 / 1 / 2 / 3 / >=4" else paste0("Arm log2 threshold +/-",parent$selected_threshold))
    result[[i]] <- make_triptych(left,middle,right,title)
    scatter <- ggplot(pdo,aes(x,y))+geom_point(size=3.5,colour="#285E79")+
      geom_text(aes(label=sub("_Untreated_PDO$","",sample)),vjust=-0.8,size=4)+
      labs(title=title,subtitle=stat_text(pdo_stat,TRUE),x="Original WES segment log2 ratio",y="Sample-mean scATLAS MP UCell activity",
        caption="Six independent untreated PDOs; existing WGS-mode CNVkit output from WES. No regression fit or RNA-guided CNA correction.")+
      theme_classic(base_size=18)+theme(plot.title=element_text(face="bold",size=19),plot.caption=element_text(size=11),
        plot.margin=margin(18,28,18,22))+scale_y_continuous(expand=expansion(mult=c(0.1,0.2)))+
      scale_x_continuous(expand=expansion(mult=c(0.2,0.2)))
    scatter_pages[[paste(kind,i)]] <<- scatter
  }
  result
}
gene_box_pages <- pair_pages("gene"); arm_box_pages <- pair_pages("arm")
write.csv(bind_rows(box_sources),file.path(tiers[["tables"]],"Auto_boxplot_source.csv"),row.names=FALSE)
write.csv(bind_rows(comparison_rows),file.path(tiers[["tables"]],"Auto_group_wilcoxon_tests.csv"),row.names=FALSE)
fixed_summary <- bind_rows(pair_summaries)
write.csv(fixed_summary,file.path(tiers[["tables"]],"Auto_fixed_pair_validation_summary.csv"),row.names=FALSE)
pdo_gene_tests <- gene_tests[gene_tests$dataset=="PDO bulk WES CNA",]
pdo_arm_tests <- arm_tests[arm_tests$dataset=="PDO bulk WES CNA",]
result_summary <- data.frame(metric=c("independent_PDO_donors","RNA_cells","scATLAS_MPs","gene_tests","arm_tests",
  "gene_FDR_below_0.05","arm_FDR_below_0.05","fixed_gene_pairs_matching_TCGA_direction","fixed_arm_pairs_matching_TCGA_direction"),
  value=c(6L,nrow(cells),17L,nrow(pdo_gene_tests),nrow(pdo_arm_tests),sum(pdo_gene_tests$fdr<0.05,na.rm=TRUE),
    sum(pdo_arm_tests$fdr<0.05,na.rm=TRUE),sum(fixed_summary$same_direction_tcga[fixed_summary$kind=="gene"],na.rm=TRUE),
    sum(fixed_summary$same_direction_tcga[fixed_summary$kind=="arm"],na.rm=TRUE)))
write.csv(result_summary,file.path(tiers[["tables"]],"Auto_result_summary.csv"),row.names=FALSE)
bundle <- list(
  list(name="Auto_gene_cna_mp_rectangles.pdf",pages=gene_pages,width=30,height=11),
  list(name="Auto_arm_cna_mp_rectangles.pdf",pages=arm_pages,width=38,height=11),
  list(name="Auto_gene_cna_mp_boxplots.pdf",pages=gene_box_pages,width=30,height=10),
  list(name="Auto_arm_cna_mp_boxplots.pdf",pages=arm_box_pages,width=30,height=10),
  list(name="Auto_PDO_cna_mp_scatterplots.pdf",pages=scatter_pages,width=14,height=9))
saveRDS(bundle,figure_cache)
saveRDS(list(manifest=manifest,scores=scores,gene_values=gene_values,arm_values=arm_values,
  gene_tests=gene_tests,arm_tests=arm_tests,mp_sets=mp_sets,reference_pairs=list(gene=reference$top_gene,arm=reference$top_arm),
  permutations=perms,threshold=parent$selected_threshold),file.path(tiers[["intermediate"]],"Auto_analysis_results.rds"))
save_figures(bundle)

# Independent numerical checks against persisted sources and a second test formulation.
check_gene <- "MYC"; check_sample <- samples[1]
g <- go[go$gene==check_gene,]; s <- segments[segments$sample==check_sample & segments$chromosome==g$chr &
  segments$start<=g$start-1 & segments$end>=g$end,]
stopifnot(nrow(s)==1L,abs(gene_values$value[gene_values$sample==check_sample & gene_values$gene==check_gene]-s$log2)<1e-12)
for(mp in mp_order) stopifnot(max(abs(scores[[mp]]-vapply(samples,function(x)mean(cells[cells$sample==x,mp]),numeric(1))))<1e-12)
for(kind in c("gene","arm")) {
  tests <- if(kind=="gene") gene_tests else arm_tests
  t <- tests[tests$dataset=="PDO bulk WES CNA",]
  stopifnot(nrow(t)==if(kind=="gene")340L else 782L,all(t$n==6L),!anyDuplicated(t[c(kind,"mp")]),
    all(t$p_value[is.finite(t$p_value)]>=1/720),all.equal(t$fdr,p.adjust(t$p_value,"BH",n=nrow(t)))==TRUE)
  vals <- if(kind=="gene") gene_values else arm_values
  for(i in seq_len(nrow(t))) {
    x <- vals[vals[[kind]]==t[[kind]][i],]; x <- x$value[match(samples,x$sample)]
    if(is.finite(t$rho[i])) stopifnot(abs(t$rho[i]-cor(x,scores[[t$mp[i]]],method="spearman"))<1e-12)
  }
}
writeLines(c("PASS: six unique donor-level DNA/RNA pairs; no treated-sample duplication.",
  "PASS: exact 17 reference MP gene sets; sample means independently reproduced.",
  "PASS: 20 genes and 46 arms; no missing CNA imputed as neutral.",
  "PASS: original segment MYC lookup reproduced independently.",
  "PASS: all PDO Spearman rho values independently checked; 720 permutations per exact two-sided P.",
  "PASS: BH correction over all 340 gene/MP and 782 arm/MP hypotheses.",
  "PASS: boxplot pairs fixed by existing reference top-ten lists, not PDO effects.",
  "LIMITATION: existing Sarek wes=false/WGS-mode CNA processing retained at user instruction.",
  "LIMITATION: relative segment CNA, estimated CN categories, small n; exploratory, no causal inference.",
  "LIMITATION: scATLAS gene predictor is RNA, whereas PDO/TCGA gene predictors are DNA CNA.",
  "LIMITATION: sex-chromosome comparisons are displayed as in reference; sex/purity/ploidy are not adjusted."),
  file.path(tiers[["reports"]],"Auto_numerical_validation.txt"))
audit_pdf()
pdo_write_run_summary(script_path,out_dir,inputs=inputs,outputs=vapply(bundle,function(x)x$name,character(1)),
  parameters=list(n_donors=6,mp_n=17,gene_n=20,arm_n=46,permutations=720,maxRank=1500,
    cna_method="original unshifted segments; no reprocessing",unit="one untreated RNA sample per DNA donor",
    caveat="WGS-mode CNVkit on WES; exploratory",pair_selection="fixed existing reference top-ten lists"))
message("Completed and verified: ",out_dir)
