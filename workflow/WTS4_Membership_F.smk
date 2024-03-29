# WTS4_Membership_F
###################################################
# normalize pattern
# normalize_pattern = ["normalize_1"]
normalize_pattern = ["n1_24sample_add3","n1_24sample_add8","n1_24sample_add25","n1_24sample_add20","n1_27sample_rm20","n1_28sample"]

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

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
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData', 
            dist=dist_data,
            range=time_range,
            N_cls=N_CLUSTERS,
            normalize_P=normalize_pattern
            )
        
rule WTS4_Membership_F:
    input:
        'output/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.RData'
    output:
        membership ='output/WTS4/{normalize_P}/{range}/{dist}/Membership_F/k_Number_{N_cls}.RData',
        cluster ='output/WTS4/{normalize_P}/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/Membership_F/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/Membership_F/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Membership_F.sh {wildcards.N_cls} {input} {output.membership} {output.cluster} >& {log}'