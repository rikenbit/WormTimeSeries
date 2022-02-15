# WTS4_Silhouette_Sample
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["4"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Eval_sample/Silhouette/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Sil_plot/Sil_gg.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            )
        
rule WTS4_Silhouette_Sample:
    input:
        sample_cls = 'output/WTS4/normalize_1/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData'
    output:
        value = 'output/WTS4/normalize_1/{range}/{dist}/Eval_sample/Silhouette/k_Number_{N_cls}.RData',
        gg_object ='output/WTS4/normalize_1/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Sil_plot/Sil_gg.RData'
    params:
        distance = 'output/WTS4/normalize_1/{range}/{dist}/Distance',
        plot = 'output/WTS4/normalize_1/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Sil_plot',
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Eval_sample/Silhouette/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/silhouette_usedist:20220215"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Eval_sample/Silhouette/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Silhouette_Sample.sh {params.distance} {input.sample_cls} {output.value} {params.plot} {output.gg_object} >& {log}'