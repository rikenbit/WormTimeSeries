# WTS4_DimReduc
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))

# Distance Data
dist_data = ["EUCL","SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]

# Dimensionality Reduction Method
DimReduc = ["tsne","umap"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{DR}_plot_docker.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc
            )
rule DimReduc_MCMI:
    input:
        m_distance = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_distance_docker.RData',
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_cls_docker.RData'
    output:
        'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{DR}_plot_docker.png'
    params:
        NL = 'data/igraph/Fig1_HNS.RData',
        EL = 'data/WTS4_Eval_behavior_fix.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{DR}_plot_docker.txt'
    container:
        "docker://yamaken37/dimreduc:20211213"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{DR}_plot_docker.log'
    shell:
        'src/WTS4_DimReduc_MCMI.sh {input.m_distance} {input.m_cls} {wildcards.DR} {output} {params.NL} {params.EL} >& {log}'