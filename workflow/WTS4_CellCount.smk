# WTS4_CellCount
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
        
rule WTS4_CellCount:
    output:
        'output/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.RData'
    params:
        dist_matrix_path = 'output/WTS4/{normalize_P}/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/Distance/CellCount.log'
    shell:
        'src/WTS4_CellCount.sh {params.dist_matrix_path} {output}>& {log}'