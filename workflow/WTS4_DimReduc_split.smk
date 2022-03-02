# WTS4_DimReduc_split
###################################################
# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["5"]

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
# ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
ReClustering_method = ["CSPA"]

# Dimensionality Reduction Method
# DimReduc = ["tsne","umap"]
DimReduc = ["tsne"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_cls.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_NT.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_eval_label.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_count_sum.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_cell_count.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc
            ),
rule WTS4_DimReduc_split:
    input:
        m_distance = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_distance/k_Number_{N_cls}.RData',
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData'
    output:
        gg_cls = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_cls.png',
        gg_NT = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_NT.png',
        gg_eval_label = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_eval_label.png',
        gg_count_sum = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_count_sum.png',
        gg_cell_count = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}_gg_cell_count.png'
    params:
        NL = 'data/igraph/Fig1_HNS.RData',
        EL = 'data/WTS4_Eval_behavior_fix.xlsx',
        cell_count = 'output/WTS4/normalize_1/{range}/{dist}/Distance/CellCount.RData',
        count_sum  = 'output/WTS4/normalize_1/{range}/{dist}/ClsCount/k_Number_{N_cls}/df_count_sum.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/dimreduc:20211213"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_{DR}/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_DimReduc_split.sh {input.m_distance} {input.m_cls} {wildcards.DR} {params.NL} {params.EL} {params.cell_count} {params.count_sum} {output.gg_cls} {output.gg_NT} {output.gg_eval_label} {output.gg_count_sum} {output.gg_cell_count} >& {log}'
