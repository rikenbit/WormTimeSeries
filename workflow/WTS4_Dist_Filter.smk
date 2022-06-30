# WTS4_Dist_Filter
###################################################
# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# normalize pattern
# normalize_pattern = ["normalize_1"]
normalize_pattern = ["n1_24sample_add3","n1_24sample_add8","n1_24sample_add25","n1_24sample_add20","n1_27sample_rm20","n1_28sample"]

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
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.txt'
    container:
        "docker://yamaken37/dist_filter:20220624"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/Distance/Ds_F.log'
    shell:
        'src/WTS4_Dist_Filter.sh {input} {output}>& {log}'