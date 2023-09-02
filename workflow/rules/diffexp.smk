rule count_matrix:
    input:
        expand(
            "results/star/{unit.sample_name}_{unit.unit_name}/ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
    output:
        "results/counts/all.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample_name"].tolist(),
        strand=get_strandedness(units),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"

rule normalize_counts:
    input:
        counts="results/counts/all.tsv",
        gtf="resources/genome.gtf",
    output:
        tpm="results/counts/normalized_tpm.tsv",
        cpm="results/counts/normalized_cpm.tsv",
    log:
        "logs/normalize_counts.log",
    script:
        "../scripts/rnanorm.py"

rule gene_2_symbol:
    input:
        counts="{prefix}.tsv",
    output:
        symbol="{prefix}.symbol.tsv",
    params:
        species=get_bioc_species_name(),
    log:
        "logs/gene2symbol/{prefix}.log",
    conda:
        "../envs/biomart.yaml"
    script:
        "../scripts/gene2symbol.R"


rule deseq2_init:
    input:
        counts="results/counts/all.tsv",
    output:
        "results/deseq2/all.rds",
        "results/deseq2/normcounts.tsv",
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "results/deseq2/all.rds",
    output:
        report("results/pca.{variable}.svg", "../report/pca.rst"),
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.{variable}.log",
    script:
        "../scripts/plot-pca.R"

rule heatmap:
    input:
        "results/deseq2/all.rds",
    output:
        report("results/heatmap.{variable}.svg", "../report/heatmap.rst"),
    params:
        heatmap_labels=config["heatmap"]["labels"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/heatmap.{variable}.log",
    script:
        "../scripts/plot-heatmap.R"

rule volcano:
    input:
        "results/deseq2/all.rds",
    output:
        report("results/diffexp/{contrast}.volcanoplot.svg", "../report/volcanoplot.rst"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/{contrast}.volcanoplot.log",
    script:
        "../scripts/plot-volcano.R"


rule deseq2:
    input:
        "results/deseq2/all.rds",
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2.R"

rule deseq2_expressiontable:
    input:
        "results/deseq2/all.rds",
    output:
        table=report("results/diffexp/{contrast}.expressiontable.tsv", "../report/expressiontable.rst"),
    params:
        contrast=get_contrast,
        samples=config["samples"],
        model=config["diffexp"]["model"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.expresssiontable.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-exptable.R"
