# WTS4_ReClustering_noisy
###################################################
# NOISE SAMPLE PATERN
NOISE_TEST = ["n1_28sample","n1_24sample_add3","n1_24sample_add8","n1_24sample_add20","n1_24sample_add25"]
# NOISE_TEST = ["n1_28sample"]

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter","all"]
# time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["MCMIHOOI"]

rule all:
    input:
        expand('output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            NOISE=NOISE_TEST
            ),
        expand('output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_distance/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            NOISE=NOISE_TEST
            ),
        expand('output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            NOISE=NOISE_TEST
            )
        
rule WTS4_ReClustering_noisy:
    input:
        Mem_matrix = 'output/WTS4/{NOISE}/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
    output:
        m_data = 'output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.RData',
        m_distance = 'output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_distance/k_Number_{N_cls}.RData',
        m_cls = 'output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData'
    benchmark:
        'benchmarks/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/reclustering:20211118"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_ReClustering.sh {wildcards.N_cls} {wildcards.Re_cls} {input.Mem_matrix} {output.m_data} {output.m_distance} {output.m_cls} >& {log}'