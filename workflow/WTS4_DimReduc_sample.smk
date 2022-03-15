# WTS4_DimReduc_sample
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["5"]

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# Dimensionality Reduction Method
# DimReduc = ["tsne","umap"]
DimReduc = ["tsne"]

# data time range
time_range = ["stimAfter"]

# normalize pattern
# normalize_pattern = ["normalize_1"]
normalize_pattern = ["n1_28sample"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/{DR}/table.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            DR=DimReduc,
            normalize_P=normalize_pattern
            )
rule WTS4_DimReduc_sample:
    input:
        m_data = 'output/WTS4/{normalize_P}/{range}/{dist}/MCMIHOOI/Merged_data/k_Number_{N_cls}.RData',
        sample_cls = 'output/WTS4/{normalize_P}/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData'
    output:
        'output/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/{DR}/table.png'
    params:
        NL = 'data/igraph/Fig1_HNS.RData',
        EL = 'data/WTS4_Eval_behavior_fix.xlsx',
        input_path = 'output/WTS4/{normalize_P}/{range}/{dist}/Distance',
        cell_count = 'output/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.RData'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/{DR}/table.txt'
    container:
        "docker://yamaken37/dimreduc_mcmi:20211224"
        # コンテナはDimReduc_MCMIのを使用
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/{DR}/table.log'
    shell:
        'src/WTS4_DimReduc_sample.sh {input.m_data} {output} {params.input_path} {wildcards.DR} {wildcards.N_cls} {params.NL} {params.EL} {input.sample_cls} {params.cell_count} >& {log}'