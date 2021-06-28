# WTS3 load
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# Distance Data
dist_data = ["EUCL","DTW","SBD"]

# Dimensionality Reduction Method
Dim_Method = ["tsne","umap"]
# Dim_Method = ["tsne"]

# Clustering Evaluation Method
cls_eval = ["purity", "ARI", "Fmeasure", "Entropy"]
# cls_eval = ["ARI"]

# data time range
time_range = ["all"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/cls_plot/SampleNumber_{N}.png', 
            N=N_SAMPLES, 
            eval=cls_eval, 
            dim_method=Dim_Method,
            dist=dist_data,
            range=time_range),
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/cls_tempdata/SampleNumber_{N}.RData', 
            N=N_SAMPLES, 
            eval=cls_eval, 
            dim_method=Dim_Method,
            dist=dist_data,
            range=time_range)

rule WTS3_Visualization:
    input:
        RData = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/{dist}.RData'
    output:
        png = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/cls_plot/SampleNumber_{N}.png',
        RData = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/cls_tempdata/SampleNumber_{N}.RData'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/cls_plot/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS3_st_mclust.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/cls_plot/SampleNumber_{N}.log'
    shell:
        'src/WTS3_vis.sh {wildcards.N} {output.png} {output.RData} {input.RData} {wildcards.eval} {wildcards.dim_method} {wildcards.range} >& {log}'
###################################################