log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("pheatmap")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[c(snakemake@wildcards[["variable"]])])
ntd <- normTransform(dds)
rownames(df) <- colnames(assay(ntd)[select,])
svg(snakemake@output[[1]])
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()
