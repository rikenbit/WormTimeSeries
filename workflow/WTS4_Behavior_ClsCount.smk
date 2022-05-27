# WTS4_Behavior_ClsCount
###################################################
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# normalize pattern
# normalize_pattern = ["normalize_1"]
normalize_pattern = ["n1_28sample"]

rule all:
    input: 
        expand('output/WTS4/{normalize_P}/{range}/{dist}/ClsCount/k_Number_{N_cls}/df_count_all.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/{range}/{dist}/ClsCount/k_Number_{N_cls}/df_count_sum.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            normalize_P=normalize_pattern
            )
        
rule WTS4_Behavior_ClsCount:
    input:
        sample_cls = 'output/WTS4/{normalize_P}/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData'
    output:
        count_all = 'output/WTS4/{normalize_P}/{range}/{dist}/ClsCount/k_Number_{N_cls}/df_count_all.RData',
        count_sum = 'output/WTS4/{normalize_P}/{range}/{dist}/ClsCount/k_Number_{N_cls}/df_count_sum.RData'
    params:
        behavior_label_path = 'data/WTS4_Eval_behavior_newlabel.xlsx',
        sample_path  = 'output/WTS4/{normalize_P}/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/ClsCount/k_Number_{N_cls}/Behavior_ClsCount.txt'
    container:
        "docker://yamaken37/eval_sample:20220106"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/ClsCount/k_Number_{N_cls}/Behavior_ClsCount.log'
    shell:
        'src/WTS4_Behavior_ClsCount.sh {input.sample_cls} {output.count_all} {output.count_sum} {wildcards.N_cls} {params.behavior_label_path} {params.sample_path} >& {log}'