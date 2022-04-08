# WTS4_Dist_List_Save
###################################################
# Distance Data
dist_data = ["EUCL","SBD_abs"]

# data time range
# time_range = ["all","stimAfter"]
time_range = ["stimAfter"]

# NORMALIZE & SAMPLES
NOISE_TEST = ["n1_28sample"]

N_SAMPLES = list(map(str, range(1, 29)))

rule all:
    input:
        expand('output/WTS4/{NOISE}/{range}/{dist}/Distance/Ds.RData', 
            dist=dist_data,
            range=time_range,
            NOISE=NOISE_TEST
            )
        
rule WTS4_Dist_List_Save:
    input:
        expand('output/WTS4/{NOISE}/{range}/{dist}/Distance/SampleNumber_{N}.RData',
            dist=dist_data,
            range=time_range,
            NOISE=NOISE_TEST,
            N=N_SAMPLES
            )
    output:
        'output/WTS4/{NOISE}/{range}/{dist}/Distance/Ds.RData'
    params:
        dist_matrix_path = 'output/WTS4/{NOISE}/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/{NOISE}/{range}/{dist}/Distance/Ds.txt'
    # conda:
    #     '../envs/myenv_WTS4_Membership.yaml'
    container:
        "docker://yamaken37/cluster_sample:20220119"
        # tidyverseが入っているコンテナ
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{NOISE}/{range}/{dist}/Distance/Ds.log'
    shell:
        'src/WTS4_Dist_List_Save.sh {params.dist_matrix_path} {output}>& {log}'