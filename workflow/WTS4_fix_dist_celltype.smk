# WTS4_fix_dist_celltype
###################################################
# N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES = ["14","26"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# data time range
time_range = ["all","stimAfter"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Distance_fix/SampleNumber_{N}.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            )
        
rule WTS4_fix_dist_celltype:
    input:
        'output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData'
    output:
        'output/WTS4/normalize_1/{range}/{dist}/Distance_fix/SampleNumber_{N}.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Distance_fix/SampleNumber_{N}.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
        # tidyverseが入っているコンテナとして選んだ
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Distance_fix/SampleNumber_{N}.log'
    shell:
        'src/WTS4_fix_dist_celltype.sh {input} {output} >& {log}'
