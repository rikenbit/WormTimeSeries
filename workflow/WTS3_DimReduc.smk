# WTS3 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
# remove not annotation ASER, not periodic ASER →17sample
N_SAMPLES.remove('4')
# N_SAMPLES.remove('7')
N_SAMPLES.remove('9')
N_SAMPLES.remove('15')
N_SAMPLES.remove('16')
N_SAMPLES.remove('22')
N_SAMPLES.remove('26')
N_SAMPLES.remove('28')

# Distance Data
dist_data = ["SBD"]

# data time range
time_range = ["all"]

# Clustering Evaluation Method
cls_eval = ["ARI"]
# cls_eval = ["purity", "ARI", "Fmeasure", "Entropy"]

# Dimensionality Reduction Method
DimReduc = ["tsne"]
# DimReduc = ["tsne","umap"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/{DR}/SampleNumber_{N}.png', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range,
            eval=cls_eval,
            DR=DimReduc
            )
        
rule DimReduc:
    input:
        RData = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/SBD.RData',
        yshift_value = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_value.RData',
        label = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/SampleNumber_{N}/label_table.RData',
        cutree = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/SampleNumber_{N}/cutree_table.RData'
    output:
        DimReduc = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{DR}/SampleNumber_{N}.png',
    params:
        args_igraph = 'data/igraph/Fig1_HNS.RData',
        args_periodic = 'output/WTS2/WTS2_PeriodicACF.csv',
        args_shift = 'ASER'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/{eval}/{DR}/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS3_DimReduc.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/{eval}/{DR}/SampleNumber_{N}.log'
    shell:
        'src/WTS3_DimReduc.sh {wildcards.N} {input.RData} {input.yshift_value} {input.label} {input.cutree} {wildcards.DR} {wildcards.eval} {output.DimReduc} >& {log}'
###################################################