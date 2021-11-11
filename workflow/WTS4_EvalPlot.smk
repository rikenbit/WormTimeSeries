# WTS4_Evaluation
###################################################
# data time range
time_range = ["stimAfter"]
# Distance Data
dist_data = ["EUCL","SBD_abs"]

# ReClustering Method
ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
# Evaluation Method
Evaluation_method = ["PseudoF","Connectivity"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/{Eval}/EvalPlot.png',
            range=time_range,
            dist=dist_data,
            Re_cls=ReClustering_method,
            Eval=Evaluation_method
            )
        
rule EvalPlot:
    output:
        EvalPlot = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/{Eval}/EvalPlot.png'
    params:
        input_path = 'output/WTS4/normalize_1/{range}/{dist}'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/{Eval}.txt'
    conda:
        '../envs/myenv_WTS4_EvalPlot.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/{Eval}.log'
    shell:
        'src/WTS4_EvalPlot.sh {params.input_path} {output.EvalPlot} {wildcards.Re_cls} {wildcards.Eval} >& {log}'