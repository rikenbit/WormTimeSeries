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
        expand('output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            )
        
rule WTS4_fix_dist_celltype:
    # input:
    #     RData = 'data/normalize_1/ReadData_{N}.RData'
    output:
        dist_matrix = 'output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS4_{dist}.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.log'
    shell:
        'src/WTS4_fix_dist_celltype.sh {output.dist_matrix}>& {log}'
