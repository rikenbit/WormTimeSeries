# WTS4_DimReduc_plot
###################################################
# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["6"]

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
# ReClustering_method = ["CSPA","MCMIHOOI"]
ReClustering_method = ["MCMIHOOI"]

# Dimensionality Reduction Method
DimReduc = ["tsne","umap"]
# DimReduc = ["tsne"]

# normalize pattern
normalize_pattern = ["normalize_1"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_cls.svg',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_NT.svg',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_eval_label.svg',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_count_sum.svg',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_cell_count.svg',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/label_table_k{N_cls}.csv',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            DR=DimReduc,
            normalize_P=normalize_pattern
            )

rule WTS4_DimReduc_plot:
    input:
        cord = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_cord/k_Number_{N_cls}.RData',
        m_cls = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData'
    output:
        gg_cls = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_cls.svg',
        gg_NT = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_NT.svg',
        gg_eval_label = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_eval_label.svg',
        gg_count_sum = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_count_sum.svg',
        gg_cell_count = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}_gg_cell_count.svg',
        csv = 'output/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/label_table_k{N_cls}.csv'
    params:
        NL = 'data/igraph/Fig1_HNS.RData',
        EL = 'data/WTS4_Eval_behavior_ACF.xlsx',
        cell_count = 'output/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.RData',
        count_sum  = 'output/WTS4/{normalize_P}/{range}/{dist}/ClsCount/k_Number_{N_cls}/df_count_sum.RData'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}.txt'
    container:
        # "docker://yamaken37/dimreduc:20211213"
        "docker://yamaken37/ggplot_svg:20230118"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/{Re_cls}/Merged_{DR}_plot/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_DimReduc_plot.sh {input.cord} {input.m_cls} {wildcards.DR} {params.NL} {params.EL} {params.cell_count} {params.count_sum} {output.gg_cls} {output.gg_NT} {output.gg_eval_label} {output.gg_count_sum} {output.gg_cell_count} {output.csv} >& {log}'
