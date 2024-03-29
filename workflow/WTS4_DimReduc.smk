# WTS4_DimReduc
###################################################
# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["5","9"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
# ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
ReClustering_method = ["CSPA", "MCMIHOOI"]

# Dimensionality Reduction Method
# DimReduc = ["tsne","umap"]
DimReduc = ["tsne"]

# normalize pattern
normalize_pattern = ["normalize_1"]
# normalize_pattern = ["n1_28sample"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}/label_table_k{N_cls}.csv',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            )
rule WTS4_DimReduc:
    input:
        m_distance = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_distance/k_Number_{N_cls}.RData',
        m_cls = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData',
        cell_count = 'output/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.RData',
        count_sum  = 'output/WTS4/{normalize_P}/{range}/{dist}/ClsCount/k_Number_{N_cls}/df_count_sum.RData'
    output:
        png = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}.png',
        csv = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}/label_table_k{N_cls}.csv'
    params:
        NL = 'data/igraph/Fig1_HNS.RData',
        EL = 'data/WTS4_Eval_behavior_fix.xlsx'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/dimreduc:20211213"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_DimReduc.sh {input.m_distance} {input.m_cls} {wildcards.DR} {output.png} {params.NL} {params.EL} {input.cell_count} {input.count_sum} {output.csv} >& {log}'


