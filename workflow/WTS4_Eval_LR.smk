# WTS4_Eval_LR
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
Evaluation_method = ["ARI","purity","Fmeasure","Entropy"]
# Evaluation_method = ["ARI"]


rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}/eval_result.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            Eval=Evaluation_method
            )
        
rule WTS4_Eval_LR:
    input:
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_cls.RData'
    output:
        eval_result = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}/eval_result.RData'

    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}.txt'
    container:
        "docker://yamaken37/eval_lr:20211129"
    # conda:
    #     '../envs/myenv_WTS4_Eval_LR.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}.log'
    shell:
        'src/WTS4_Eval_LR.sh {input.m_cls} {output.eval_result} {wildcards.Re_cls} {wildcards.Eval} >& {log}'