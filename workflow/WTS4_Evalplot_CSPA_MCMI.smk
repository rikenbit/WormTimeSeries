# WTS4_Evalplot_CSPA_MCMI
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
# ReClustering_method = ["CSPA"]

# Evaluation Method
Evaluation_method = ["PseudoF_nolabel", "Connectivity_nolabel", "kNN_nolabel", "ARI_behavior", "purity_behavior", "Fmeasure_behavior", "Entropy_behavior"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            )
        
rule WTS4_Evalplot_CSPA_MCMI:
    input:
        Mem_matrix = 'output/WTS4/normalize_1/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
    output:
        m_data = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.txt'
    container:
        "docker://docker_images"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Evalplot_CSPA_MCMI.sh {wildcards.Re_cls} {input.Mem_matrix} {output.m_data}>& {log}'