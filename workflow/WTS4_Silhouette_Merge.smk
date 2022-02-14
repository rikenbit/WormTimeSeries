# WTS4_Silhouette_Merge
###################################################
# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
# ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
ReClustering_method = ["CSPA"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/Silhouette/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Sil_plot/k_Number_{N_cls}.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Sil_gg/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            )
        
rule WTS4_Silhouette_Merge:
    input:
        m_distance = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_distance/k_Number_{N_cls}.RData',
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData'
    output:
        value = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/Silhouette/k_Number_{N_cls}.RData',
        plot = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Sil_plot/k_Number_{N_cls}.png',
        gg_object ='output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Sil_gg/k_Number_{N_cls}.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/Silhouette/k_Number_{N_cls}.RData.txt'
    container:
        "docker://yamaken37/silhouette:20220210"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/Silhouette/k_Number_{N_cls}.RData.log'
    shell:
        'src/WTS4_Silhouette_Merge.sh {input.m_distance} {input.m_cls} {output.value} {output.plot} {output.gg_object} >& {log}'