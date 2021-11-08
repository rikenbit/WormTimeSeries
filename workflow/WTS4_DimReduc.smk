# WTS4_DimReduc
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["4"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
# ReClustering_method = ["MCMIHOOI"]

# Dimensionality Reduction Method
DimReduc = ["tsne","umap"]
# DimReduc = ["tsne"]


rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{DR}_plot.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc
            )
        
rule DimReduc:
    input:
        m_distance = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_distance.RData',
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_cls.RData'
    output:
        'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{DR}_plot.png'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{DR}_plot.txt'
    conda:
        '../envs/myenv_WTS3_DimReduc.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{DR}_plot.log'
    shell:
        'src/WTS4_DimReduc.sh {input.m_distance} {input.m_cls} {wildcards.DR} {output} >& {log}'