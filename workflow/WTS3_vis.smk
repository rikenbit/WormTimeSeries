# WTS3 load
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# Distance Data
dist_data = ["SBD"]

# Dimensionality Reduction Method
DR_Method = ["umap"]

# option perplexity
# op_per = ["15"]
# option max_iter
# op_max_iter = ["3000", "5000"]

#Clustering Evaluation Method
cls_eval = ["purity", "ARI", "Fmeasure", "Entropy"]

rule all:
    input:
        expand('output/WTS3/{dist}/normalize_1/all/{dr_method}/{eval}/SampleNumber_{N}.png', 
            N=N_SAMPLES, 
            eval=cls_eval, 
            dr_method=DR_Method,
            dist=dist_data)

rule WTS3_Visualization:
    input:
        RData = 'output/WTS3/{dist}/normalize_1/all/SampleNumber_{N}/{dist}.RData'
    output:
        png = 'output/WTS3/{dist}/normalize_1/all/{dr_method}/{eval}/SampleNumber_{N}.png'
    benchmark:
        'benchmarks/WTS3/{dist}/normalize_1/all/{dr_method}/{eval}/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/{dist}/normalize_1/all/{dr_method}/{eval}/SampleNumber_{N}.log'
    shell:
        'src/WTS3_vis.sh {wildcards.N} {output.png} {input.RData} {wildcards.eval} {wildcards.dr_method} >& {log}'
###################################################