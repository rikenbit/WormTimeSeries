# WTS4_Eval_sample_behavior
###################################################

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# Evaluation Method
Evaluation_method = ["ARI","purity","Fmeasure","Entropy","NMI"]

rule all:
    input: 
        expand('output/WTS4/normalize_1/{range}/{dist}/Eval_sample/{Eval}_behavior/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Eval=Evaluation_method
            )
        
rule WTS4_Eval_sample_behavior:
    input:
         sample_cls =  'output/WTS4/normalize_1/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData'
    output:
        eval_result = 'output/WTS4/normalize_1/{range}/{dist}/Eval_sample/{Eval}_behavior/k_Number_{N_cls}.RData'
    params:
        behavior_label_path = 'data/WTS4_Eval_behavior_fix.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Eval_sample/{Eval}_behavior/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/eval_behavior:20220208"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Eval_sample/{Eval}_behavior/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Eval_sample_behavior.sh {input.sample_cls} {output.eval_result} {wildcards.Eval} {params.behavior_label_path} >& {log}'