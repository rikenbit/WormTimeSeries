# WTS4_Eval_sample
###################################################
# data time range
time_range = ["stimAfter"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_sample.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS
            )
        
rule WTS4_Eval_sample:
    input:
        sample_cls = 'output/WTS4/normalize_1/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData',
        merged_cls = 'output/WTS4/normalize_1/{range}/{dist}/MCMIHOOI/Merged_cls/k_Number_{N_cls}.RData',
        merged_data = 'output/WTS4/normalize_1/{range}/{dist}/MCMIHOOI/Merged_data/k_Number_{N_cls}.RData'
    output:
        'output/WTS4/normalize_1/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_sample.png'
    params:
        merged_label_path = 'data/WTS4_Eval_sample_fix.xlsx',
        input_path = 'output/WTS4/normalize_1/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_sample.txt'
    container:
        "docker://yamaken37/eval_sample:20220106"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_sample.log'
    shell:
        'src/WTS4_Eval_sample.sh {input.sample_cls} {output} {input.merged_cls} {input.merged_data} {params.input_path} >& {log}'