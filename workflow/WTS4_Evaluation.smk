# WTS4_Evaluation
###################################################
# data time range
time_range = ["stimAfter"]
# Distance Data
dist_data = ["EUCL","SBD_abs"]
# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))

# ReClustering Method
ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
# Evaluation Method
Evaluation_method = ["PseudoF","Connectivity"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}/eval_result_docker.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method,
            Eval=Evaluation_method
            )
        
rule Evaluation_docker:
    input:
        m_data = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_data_docker.RData',
        m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_cls_docker.RData'
    output:
        eval_result = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}/eval_result_docker.RData'

    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}_docker.txt'
    conda:
        '../envs/myenv_WTS4_Evaluation.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}_docker.log'
    shell:
        'src/WTS4_Evaluation.sh {input.m_data} {input.m_cls} {output.eval_result} {wildcards.Re_cls} {wildcards.Eval} >& {log}'
# rule Evaluation:
#     input:
#         m_data = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_data.RData',
#         m_cls = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/merged_cls.RData'
#     output:
#         eval_result = 'output/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}/eval_result.RData'

#     benchmark:
#         'benchmarks/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}.txt'
#     conda:
#         '../envs/myenv_WTS4_Evaluation.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/normalize_1/{range}/{dist}/{N_cls}_Clusters/{Re_cls}/{Eval}.log'
#     shell:
#         'src/WTS4_Evaluation.sh {input.m_data} {input.m_cls} {output.eval_result} {wildcards.Re_cls} {wildcards.Eval} >& {log}'