rule align:
    input:
        unpack(get_fq),
        index="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln=temp("results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam"),
        reads_per_gene="results/star/{sample}_{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}_{unit}.log",
    params:
        idx=lambda wc, input: input.index,
        extra=lambda wc, input: f'--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {input.gtf} {config["params"]["star"]}',
    threads: 24
    wrapper:
        "v1.21.4/bio/star/align"

rule samtools_index:
    input:
        "results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam.bai"
    params:
        threads=24
    shell:
        "samtools index -@ {params.threads} {input}"