# WTS4_Membership_Label
###################################################
# SAMPLE NUMBER
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
# N_SAMPLES = ["10"]

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# Distance Data
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# data group
NOISE_TEST = ["n1_noise"]

N_TRY = list(map(str, range(1, 11)))
# N_TRY = ["1"]

rule all:
    input:
        expand('output/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
            dist=dist_data,
            range=time_range,
            N=N_SAMPLES,
            NOISE=NOISE_TEST,
            TRY=N_TRY,
            N_cls=N_CLUSTERS
            )
        
rule WTS4_Membership_Label:
    input:
        Noise_matrix = 'output/WTS4/{NOISE}_{TRY}/{range}/{dist}/Distance/SampleNumber_{N}.RData'
    output:
        Mem_matrix = 'output/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
    params:
        dist_matrix_path = 'output/WTS4/normalize_1/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
        # tidyverseが入っているコンテナ
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/Membership/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Membership_Label.sh {wildcards.N_cls} {params.dist_matrix_path} {input.Noise_matrix} {output.Mem_matrix} {wildcards.N} >& {log}'