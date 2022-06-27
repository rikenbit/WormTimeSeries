# WTS4_Membership_F_F
###################################################
# normalize pattern
normalize_pattern = ["normalize_1"]

# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["3"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# data time range
# time_range = ["all","stimAfter"]
time_range = ["stimAfter"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/Membership_F/k_Number_{N_cls}.RData', 
            dist=dist_data,
            range=time_range,
            N_cls=N_CLUSTERS,
            normalize_P=normalize_pattern
            )
        
rule WTS4_Membership_F:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.RData',
            dist=dist_data,
            range=time_range,
            normalize_P=normalize_pattern
            )
    output:
        'output/WTS4/{normalize_P}/{range}/{dist}/Membership_F/k_Number_{N_cls}.RData'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/Membership_F/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/Membership_F/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Membership_F.sh {wildcards.N_cls} {input} {output}>& {log}'