# WTS3 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# Distance Data
dist_data = ["SBD"]

# Dimensionality Reduction Method
Dim_Method = ["tsne"]
# Dim_Method = ["tsne","umap"]

# Clustering Evaluation Method
cls_eval = ["ARI"]
# cls_eval = ["purity", "ARI", "Fmeasure", "Entropy"]

# filtering data
df_filter = ["stim_cell","stim_cluster"]
# data time range
time_range = ["all"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_plot/{df_f}/SampleNumber_{N}.png', 
            N=N_SAMPLES, 
            eval=cls_eval, 
            dim_method=Dim_Method,
            dist=dist_data,
            df_f=df_filter,
            range=time_range),
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_table/{df_f}/SampleNumber_{N}.RData', 
            N=N_SAMPLES, 
            eval=cls_eval, 
            dim_method=Dim_Method,
            dist=dist_data,
            df_f=df_filter,
            range=time_range)

rule WTS3_plot:
    input:
        Neuron = 'data/normalize_1/ReadData_{N}.RData',
        stim = 'data/stimulation/stim_{N}.RData',
        mCherry = 'data/mCherry/mCherry_{N}.RData',
        Position = 'data/Position/Position_{N}.RData',
        tempdata = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/cls_tempdata/SampleNumber_{N}.RData'
        
    output:
        png = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_plot/{df_f}/SampleNumber_{N}.png',
        RData = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_table/{df_f}/SampleNumber_{N}.RData'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_plot/{df_f}/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS3_plot.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_plot/{df_f}/SampleNumber_{N}.log'
    shell:
        'src/WTS3_plot.sh {wildcards.N} {input.Neuron} {input.stim} {input.mCherry} {input.Position} {input.tempdata} {output.png} {wildcards.eval} {wildcards.dim_method} {wildcards.df_f} {wildcards.range} {output.RData} >& {log}'
###################################################