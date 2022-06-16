# WTS4_Eval_Sample_Fig4

###################################################
# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["9"]
N_CLUSTERS = ["6"]

# Distance Data
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# normalize pattern
normalize_pattern = ["normalize_1"]
# normalize_pattern = ["n1_28sample"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_Sample_Fig4_k.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            normalize_P=normalize_pattern
            )
        
rule WTS4_Eval_Sample_Fig4:
    input:
        sample_cls = 'output/WTS4/{normalize_P}/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData'
    output:
        No_of_cells = 'output/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_Sample_Fig4_k.png'
    params:
        weight_path = 'output/WTS4/{normalize_P}/{range}/{dist}/MCMIHOOI/Merged_data',
        sample_path = 'output/WTS4/{normalize_P}/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_Sample_Fig4_k.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/{range}/{dist}/DimReduc_sample/k_Number_{N_cls}/Eval_Sample_Fig4_k.log'
    shell:
        'src/WTS4_Eval_Sample_Fig4.sh {params.weight_path} {params.sample_path} {input.sample_cls} {output.No_of_cells} {wildcards.N_cls} >& {log}'