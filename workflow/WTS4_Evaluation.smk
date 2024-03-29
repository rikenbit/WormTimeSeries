# WTS4_Evaluation
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
# ReClustering_method = ["CSPA"]

# Evaluation Method
Evaluation_method = ["PseudoF","Connectivity","kNN"]
# Evaluation_method = ["KNN"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/{Eval}/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            Eval=Evaluation_method
            )
        
rule Evaluation:
    input:
        m_data = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data/k_Number_{N_cls}.RData',
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_cls/k_Number_{N_cls}.RData'
    output:
        eval_result = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/{Eval}/k_Number_{N_cls}.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/{Eval}/k_Number_{N_cls}.txt'
    # conda:
    #     '../envs/myenv_WTS4_Evaluation.yaml'
    container:
        "docker://yamaken37/evaluation:20220117"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval/{Eval}/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Evaluation.sh {input.m_data} {input.m_cls} {output.eval_result} {wildcards.Re_cls} {wildcards.Eval} {wildcards.N_cls} >& {log}'