# WTS4_yshift_merge
###################################################
# Distance Data
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# input matrix
input_matrix = ["ave"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Shift_FM_{in_mat}/SampleNumber_ALL.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range,
            in_mat=input_matrix
            )
        
        
rule WTS4_yshift_merge:
    output:
        Shift = 'output/WTS4/normalize_1/{range}/{dist}/Shift_FM_{in_mat}/SampleNumber_ALL.RData'
    params:
        mat_dir = 'output/WTS4/normalize_1/{range}/{dist}/Shift_F'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Shift_FM_{in_mat}/SampleNumber_ALL.txt'
    container:
        "docker://docker_images"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Shift_FM_{in_mat}/SampleNumber_ALL.log'
    shell:
        'src/WTS4_yshift_merge.sh {params.mat_dir} {output.Shift}>& {log}'