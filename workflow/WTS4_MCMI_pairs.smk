# WTS4_MCMI_pairs
###################################################
# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["9"]
# N_CLUSTERS = list(map(str, range(2, 14)))

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["MCMIHOOI"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Cluster.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Neuron_type.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Class.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            )
        
rule WTS4_MCMI_pairs:
    input:
        Merged_data = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.RData',
        Consistency = 'output/WTS4/normalize_1/{range}/{dist}/ClsCount/k_Number_{N_cls}/df_count_sum.RData',
        No_of_cells = 'output/WTS4/normalize_1/{range}/{dist}/Distance/CellCount.RData',
        Cluster = 'output/WTS4/normalize_1/{range}/{dist}/MCMIHOOI/Merged_cls/k_Number_{N_cls}.RData'

    output:
        Cluster = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Cluster.png',
        Neuron_type = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Neuron_type.png',
        Class = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Class.png'
    params:
        NL = 'data/igraph/Fig1_HNS.RData',
        EL = 'data/WTS4_Eval_behavior_fix.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/MCMI_pairs.txt'
    container:
        "docker://yamaken37/mcmi_pairs:20220309"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/MCMI_pairs.log'
    shell:
        'src/WTS4_MCMI_pairs.sh {input.Merged_data} {input.Consistency} {input.No_of_cells} {input.Cluster}  {params.NL} {params.EL} {output.Cluster} {output.Neuron_type} {output.Class}>& {log}'