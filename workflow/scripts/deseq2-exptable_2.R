log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

norm_counts = counts(dds, normalized=T)

coldata <- read.table(snakemake@params[["samples"]], header=TRUE, check.names=FALSE)
coldata <- coldata[order(row.names(coldata)), , drop=F]
condition <- snakemake@params[["model"]]
condition <- sub("~", "", condition)
res_tidy <- results(dds, tidy = TRUE)
for (i in unique(coldata[[condition]]))
{
    mean <- data.frame(mean=rowMeans(norm_counts[,c(coldata[coldata[[condition]]==i,]$sample_name)]))
    mean <- cbind(rownames(mean), mean)
    rownames(mean) <- NULL
    colnames(mean) <- c("row", paste("baseMean_", i, sep=""))
    res_tidy<-merge(res_tidy, mean)
}

norm_counts <- data.frame(row = row.names(norm_counts), norm_counts, row.names = NULL)

# sort by p-value
res_tidy <- res_tidy[order(res_tidy$padj),]
res_tidy<-merge(res_tidy, norm_counts, by='row')

# store results
colnames(res_tidy)[1] = "gene"
write.table(res_tidy, file=snakemake@output[["table"]], row.names=FALSE, sep='\t')
