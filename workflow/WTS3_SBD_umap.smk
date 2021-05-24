# WTS3 SBD load
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# Dimensionality Reduction Method
DR_Method = ["umap"]

# option perplexity
# op_per = ["15"]
# option max_iter
# op_max_iter = ["3000", "5000"]

#clustering evaluation method
cls_eval = ["purity"]

rule all:
    input:
        expand('output/WTS3/SBD/normalize_1/all/{dr_method}/{eval}/SampleNumber_{N}.png', N=N_SAMPLES, eval=cls_eval, dr_method=DR_Method)

rule SBD_Purity:
    input:
        RData = 'output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.RData'
    output:
        png = 'output/WTS3/SBD/normalize_1/all/{dr_method}/{eval}/SampleNumber_{N}.png'
    benchmark:
        'benchmarks/WTS3/SBD/normalize_1/all/{dr_method}/{eval}/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/SBD/normalize_1/all/{dr_method}/{eval}/SampleNumber_{N}.log'
    shell:
        'src/WTS3_SBD_umap.sh {wildcards.N} {output.png} {input.RData} {wildcards.eval} {wildcards.dr_method} >& {log}'
###################################################