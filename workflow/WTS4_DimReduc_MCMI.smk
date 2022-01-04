# WTS4_DimReduc
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["5"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# Dimensionality Reduction Method
DimReduc = ["tsne","umap"]
# DimReduc = ["tsne"]

# data time range
time_range = ["stimAfter"]
# ReClustering method
ReClustering_method = ["MCMIHOOI"]

rule all:
    input:
        # output/WTS4/normalize_1/stimAfter/SBD_abs/5_Clusters/DimReduc_MCMI/tsne/table.png
        expand('output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/DimReduc_{Re_cls}/{DR}/table.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc
            )
rule DimReduc_MCMI:
    input:
        #output/WTS4/normalize_1/stimAfter/SBD_abs/5_Clusters/MCMIHOOI/merged_data.RData
        m_data = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_data.RData'
    output:
        'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/DimReduc_{Re_cls}/{DR}/table.png'
    params:
        NL = 'data/igraph/Fig1_HNS.RData',
        EL = 'data/WTS4_Eval_behavior_fix.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/DimReduc_{Re_cls}/{DR}/table.txt'
    container:
        "docker://yamaken37/dimreduc_mcmi:20211224"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/DimReduc_{Re_cls}/{DR}/table.log'
    shell:
        'src/WTS4_DimReduc_MCMI.sh {input.m_data} {output} {wildcards.dist} {wildcards.DR} {wildcards.N_cls} {params.NL} {params.EL} >& {log}'