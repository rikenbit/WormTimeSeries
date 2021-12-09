# WTS4_EvalPlot
###################################################
# data time range
time_range = ["stimAfter"]
# Distance Data
dist_data = ["EUCL","SBD_abs"]

# ReClustering Method
ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]
# Evaluation Method
# Evaluation_method = ["PseudoF","Connectivity"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_behavior.png',
            range=time_range,
            dist=dist_data,
            Re_cls=ReClustering_method
            )
rule EvalPlot_behavior:
    output:
        EvalPlot = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_behavior.png'
    params:
        input_path = 'output/WTS4/normalize_1/{range}/{dist}'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_behavior.txt'
    conda:
        '../envs/myenv_WTS4_EvalPlot.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_behavior.log'
    shell:
        'src/WTS4_EvalPlot.sh {params.input_path} {output.EvalPlot} {wildcards.Re_cls} >& {log}'

# rule EvalPlot_docker:
#     output:
#         EvalPlot = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_docker.png'
#     params:
#         input_path = 'output/WTS4/normalize_1/{range}/{dist}'
#     benchmark:
#         'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_docker.txt'
#     conda:
#         '../envs/myenv_WTS4_EvalPlot.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot_docker.log'
#     shell:
#         'src/WTS4_EvalPlot.sh {params.input_path} {output.EvalPlot} {wildcards.Re_cls} >& {log}'

# rule EvalPlot:
#     output:
#         EvalPlot = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot.png'
#     params:
#         input_path = 'output/WTS4/normalize_1/{range}/{dist}'
#     benchmark:
#         'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot.txt'
#     conda:
#         '../envs/myenv_WTS4_EvalPlot.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/EvalPlot.log'
#     shell:
#         'src/WTS4_EvalPlot.sh {params.input_path} {output.EvalPlot} {wildcards.Re_cls} >& {log}'