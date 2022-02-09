# WTS4_Eval_behavior
###################################################
# data time range
time_range = ["stimAfter"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# ReClustering Method
ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
# ReClustering_method = ["MCMIHOOI"]

# Evaluation Method
Evaluation_method = ["ARI_behavior","purity_behavior","Fmeasure_behavior","Entropy_behavior","NMI_behavior","AMI_behavior"]
# Evaluation_method = ["NMI_behavior"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/{Eval}/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            Eval=Evaluation_method
            )
        
rule WTS4_Eval_behavior:
    input:
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData'
    output:
        eval_result = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/{Eval}/k_Number_{N_cls}.RData'
    params:
        behavior_label_path = 'data/WTS4_Eval_behavior_fix.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/{Eval}/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/eval_behavior:20220208"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/{Eval}/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Eval_behavior.sh {input.m_cls} {output.eval_result} {wildcards.Re_cls} {wildcards.Eval} {params.behavior_label_path} >& {log}'