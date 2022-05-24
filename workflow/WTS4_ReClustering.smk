# WTS4_ReClustering 
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
# ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
ReClustering_method = ["CSPA","MCMIHOOI"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_distance/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            )
        
rule ReClustering:
    input:
        Mem_matrix = 'output/WTS4/normalize_1/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
    output:
        m_data = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.RData',
        m_distance = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_distance/k_Number_{N_cls}.RData',
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/reclustering:20211118"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_ReClustering.sh {wildcards.N_cls} {wildcards.Re_cls} {input.Mem_matrix} {output.m_data} {output.m_distance} {output.m_cls} >& {log}'