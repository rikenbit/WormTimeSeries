# WTS4_DimReduc_cord
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["9"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["CSPA","MCMIHOOI"]
# ReClustering_method = ["MCMIHOOI"]

# Dimensionality Reduction Method
DimReduc = ["tsne","umap"]
# DimReduc = ["tsne"]

# normalize pattern
normalize_pattern = ["normalize_1"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_cord/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            )
rule WTS4_DimReduc_cord:
    input:
        'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_distance/k_Number_{N_cls}.RData'
    output:
        'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_cord/k_Number_{N_cls}.RData'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_cord/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/dimreduc:20211213"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_cord/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_DimReduc_cord.sh {input} {wildcards.DR} {output} >& {log}'
