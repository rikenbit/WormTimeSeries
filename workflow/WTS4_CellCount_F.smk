# WTS4_CellCount_F
###################################################
# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# normalize pattern
normalize_pattern = ["normalize_1"]
# normalize_pattern = ["n1_28sample"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.RData',
            range=time_range,
            dist=dist_data,
            normalize_P=normalize_pattern
            )
        
rule WTS4_CellCount_F:
    input:
        'output/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.RData'
    output:
        'output/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.RData'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.log'
    shell:
        'src/WTS4_CellCount_F.sh {input} {output}>& {log}'