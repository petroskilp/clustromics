rule parse_counts:
    input:
        normcounts="results/counts/normalized_tpm.tsv",
    output:
        parsedcounts="results/ml/parsed_tpm.tsv",
    log:
        "logs/ml/parse_counts.log",
    script:
        "../scripts/parsecounts.py"

rule importance:
    input:
        parsedcounts="results/ml/parsed_tpm.tsv",
    output:
        classes="results/ml/{variable}_classes.npy",
        importance="results/ml/{variable}_importance.csv",
    log:
        "logs/ml/importance_{variable}.log",
    conda:
        "../envs/ml.yaml"
    script:
        "../scripts/importance.py"

rule ml:
    input:
        parsedcounts="results/ml/parsed_tpm.tsv",
        classes="results/ml/{variable}_classes.npy",
    output:
        model="results/ml/{variable}_{model_name}_model.pkl",
        roc="results/ml/{variable}_{model_name}_roc_curve.png",
    log:
        "logs/ml/models_{variable}_{model_name}.log",
    conda:
        "../envs/ml.yaml"
    script:
        "../scripts/ml.py"