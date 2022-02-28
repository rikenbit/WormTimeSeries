# WTS4_noisy
###################################################

# WTS4_Membership
###################################################
# NOISE SAMPLE PATERN
NOISE_TEST = ["n1_28sample","n1_24sample_add3","n1_24sample_add8","n1_24sample_add20","n1_24sample_add25"]

# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = list(map(str, range(2, 4)))

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]
# data time range
# time_range = ["all","stimAfter"]
time_range = ["stimAfter"]


rule all:
    input:
        expand('output/WTS4/{NOISE_TEST}/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
            dist=dist_data,
            range=time_range,
            N_cls=N_CLUSTERS
            )
        
rule WTS4_Membership:
    input:
        expand('output/WTS4/{NOISE_TEST}/{range}/{dist}/Distance/SampleNumber_{N}.RData',
            dist=dist_data,
            range=time_range,
            N=N_SAMPLES
            )
    output:
        Mem_matrix = 'output/WTS4/{NOISE_TEST}/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
    params:
        dist_matrix_path = 'output/WTS4/{NOISE_TEST}/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/{NOISE_TEST}/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
    # conda:
    #     '../envs/myenv_WTS4_Membership.yaml'
    container:
        "docker://yamaken37/cluster_sample:20220119"
        # tidyverseが入っているコンテナ
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{NOISE_TEST}/{range}/{dist}/Membership/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Membership.sh {wildcards.N_cls} {params.dist_matrix_path} {output.Mem_matrix}>& {log}'