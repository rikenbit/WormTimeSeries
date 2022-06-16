# WTS4_dhyper
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["5","9"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["CSPA", "MCMIHOOI"]

# Dimensionality Reduction Method
DimReduc = ["tsne"]

# normalize pattern
normalize_pattern = ["normalize_1"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_dhyper/p_table_k{N_cls}.csv',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_dhyper/q_table_k{N_cls}.csv',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            )
rule WTS4_dhyper:
    input:
        csv = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/label_table_k{N_cls}.csv'
    output:
        csv = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_dhyper/p_table_k{N_cls}.csv',
        csv_q = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_dhyper/q_table_k{N_cls}.csv'
    params:
        ann = 'data/label_table_ann.csv'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_dhyper/p_table_k{N_cls}.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_dhyper/p_table_k{N_cls}.log'
    shell:
        'src/WTS4_dhyper.sh {input.csv} {params.ann} {output.csv} {output.csv_q} >& {log}'
