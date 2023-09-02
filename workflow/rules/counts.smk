rule feature_counts:
    input:
        bam="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bai="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam.bai",
    output:
        counts="results/feature_counts/{sample}_{unit}.txt",
        summary="results/qc/feature_counts/{sample}_{unit}.summary"
    params:
        annotation="resources/genome.gtf",
        extra=""#feature_counts_params(sample)
    threads:
        24
    log:
        "logs/feature_counts/{sample}_{unit}.txt"
    shell:
        "featureCounts "
        "{params.extra} "
        "-a {params.annotation} "
        "-o {output.counts} "
        "-T {threads} "
        "{input.bam} > {log} 2>&1; "
        "mv \"{output.counts}.summary\" {output.summary}"


rule merge_counts:
    input:
        [f"results/feature_counts/{sample}_{unit}.txt" for sample, unit in zip(samples["sample_name"], units["unit_name"])]

    output:
        complete="results/feature_counts/merged.tsv",
        raw="results/feature_counts/merged_raw.tsv",
        lengths="results/feature_counts/lengths.tsv"
    run:
        indices = ["Geneid", "Chr", "Start", "End", "Strand","Length"]
        frames = [
            pd.read_csv(fp, sep="\t", skiprows=1, index_col=indices)
            for fp in input
        ]
        merged = pd.concat(frames, axis=1)

        merged = merged.rename(columns=lambda c: Path(c).stem)
        merged.to_csv(output.complete, sep="\t", index=True)

        # Save counts compatible with rnanorm
        merged = merged.reset_index(level="Geneid")
        merged = merged.rename(columns={"Geneid":"FEATURE_ID"})
        merged.to_csv(output.raw, sep="\t", index=False)

        lengths = merged.reset_index(level="Length")
        lengths[["FEATURE_ID", "Length"]].to_csv(output.lengths, sep="\t", index=False)

rule normalize_counts:
    input:
        counts="results/feature_counts/merged_raw.tsv",
        lengths="results/feature_counts/lengths.tsv"
    output:
        tpm="results/feature_counts/normalized_tpm.tsv",
        cpm="results/feature_counts/normalized_cpm.tsv",
    log:
        "logs/normalize_counts.log",
    script:
        "../scripts/rnanorm.py"