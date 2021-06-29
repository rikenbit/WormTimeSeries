# WTS3 SBD 
###################################################
# N_SAMPLES = list(map(str, range(1, 29)))
# N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')

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
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/table/{df_f}/table_neuron.csv',
            eval=cls_eval, 
            dim_method=Dim_Method,
            dist=dist_data,
            df_f=df_filter,
            range=time_range),
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/table/{df_f}/table_rmSensory.csv',
            eval=cls_eval, 
            dim_method=Dim_Method,
            dist=dist_data,
            df_f=df_filter,
            range=time_range),
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_table/{df_f}/shift_table.csv',
            eval=cls_eval, 
            dim_method=Dim_Method,
            dist=dist_data,
            df_f=df_filter,
            range=time_range)

rule WTS3_table:
    output:
        csv_1 = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/table/{df_f}/table_neuron.csv',
        csv_2 = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/table/{df_f}/table_rmSensory.csv',
        shift_csv_1 = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_table/{df_f}/shift_table.csv'
    params:
        tempdata_dir = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/cls_tempdata/SampleNumber_',
        shift_dir = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/shift_table/{df_f}/SampleNumber_'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/table/{df_f}/table.txt'
    conda:
        '../envs/myenv_WTS3_plot.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/{eval}/{dim_method}/table/{df_f}/table.log'
    shell:
        'src/WTS3_table.sh {params.tempdata_dir} {wildcards.df_f} {output.csv_1} {output.csv_2}  {params.shift_dir} {output.shift_csv_1}>& {log}'
###################################################