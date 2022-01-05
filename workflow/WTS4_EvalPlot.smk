# WTS4_EvalPlot
###################################################
# data time range
time_range = ["stimAfter"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# ReClustering Method
ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
# ReClustering_method = ["CSPA"]

# label
label = ["LR","behavior"]
# label = ["LR"]
rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_{lab}.png',
            range=time_range,
            dist=dist_data,
            Re_cls=ReClustering_method,
            lab=label
            )

rule EvalPlot:
    output:
        EvalPlot = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_{lab}.png'
    params:
        input_path = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Eval'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_{lab}.txt'
    conda:
        '../envs/myenv_WTS4_EvalPlot.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_{lab}.log'
    shell:
        'src/WTS4_EvalPlot.sh {params.input_path} {output.EvalPlot} {wildcards.lab} >& {log}'