# WTS4_Dist_Filter
###################################################
# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# normalize pattern
normalize_pattern = ["normalize_1"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.RData',
            range=time_range,
            dist=dist_data,
            normalize_P=normalize_pattern
            )
        
rule WTS4_Dist_Filter:
    input:
        'output/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds.RData'
    output:
        'output/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.RData'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.RData.txt'
    container:
        "docker://yamaken37/silhouette_usedist:20220215"
        
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.RData.log'
    shell:
        'src/WTS4_Dist_Filter.sh {input} {output}>& {log}'