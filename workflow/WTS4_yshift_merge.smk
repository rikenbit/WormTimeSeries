# WTS4_yshift_merge
###################################################
# N_SAMPLES = list(map(str, range(1, 29)))

# # remove artifact
# N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')

# Distance Data
dist_data = ["SBD_abs"]
# data time range
time_range = ["stimAfter"]
# input matrix
input_matrix = ["Shift_F"]
# stat methods
stat_value = ["mean","sd","count"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{in_mat}M_{stat}/SampleNumber_ALL.RData', 
            dist=dist_data,
            range=time_range,
            in_mat=input_matrix,
            stat=stat_value
            )
        
rule WTS4_yshift_merge:
    output:
        'output/WTS4/normalize_1/{range}/{dist}/{in_mat}M_{stat}/SampleNumber_ALL.RData'
    params:
        Shift_mat = 'output/WTS4/normalize_1/stimAfter/SBD_abs/{in_mat}'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{in_mat}M_{stat}/SampleNumber_ALL.txt'
    container:
        "docker://yamaken37/yshift_visualize:20220921"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{in_mat}M_{stat}/SampleNumber_ALL.log'
    shell:
        'src/WTS4_yshift_merge.sh {params.Shift_mat} {output} {wildcards.stat} >& {log}'
