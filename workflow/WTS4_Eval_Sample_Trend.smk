# WTS4_Eval_sample
###################################################
# data time range
time_range = ["stimAfter"]

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["9"]

# normalize pattern
# normalize_pattern = ["normalize_1"]
normalize_pattern = ["n1_28sample"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_Sample_trend.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            normalize_P=normalize_pattern
            )
        
rule WTS4_Eval_Sample_Trend:
    input:
        sample_cls = 'output/WTS4/{normalize_P}/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData',
        merged_cls = 'output/WTS4/{normalize_P}/{range}/{dist}/MCMIHOOI/Merged_cls/k_Number_{N_cls}.RData',
        merged_data = 'output/WTS4/{normalize_P}/{range}/{dist}/MCMIHOOI/Merged_data/k_Number_{N_cls}.RData'
    output:
        'output/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_Sample_trend.png'
    params:
        merged_label_path = 'data/WTS4_Eval_sample_fix.xlsx',
        input_path = 'output/WTS4/{normalize_P}/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_Sample_trend.txt'
    container:
        "docker://yamaken37/eval_sample:20220106"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_Sample_trend.log'
    shell:
        'src/WTS4_Eval_Sample_Trend.sh {input.sample_cls} {output} {input.merged_cls} {input.merged_data} {params.input_path} >& {log}'