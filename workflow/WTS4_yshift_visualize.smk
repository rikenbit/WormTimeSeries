# WTS4_yshift_visualize
###################################################
N_SAMPLES = list(map(str, range(1, 29)))

# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')


# Distance Data
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# input matrix
# input_mat = ["Shift_F","Shift_F_M_ave","Shift_F_M_sd"]
input_matrix = ["Shift_F"]

# value type
value_type =["zahlen","abs"]

# filter (label combination)
label_comb =["No_F","ALL","NaCl","1n","1p","1np"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{in_mat}/SampleNumber_{N}/{v_type}/{l_comb}.png', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range,
            in_mat=input_matrix,
            v_type=value_type,
            l_comb=label_comb
            )
        
rule WTS4_yshift_visualize:
    input:
        RData = 'output/WTS4/normalize_1/{range}/{dist}/{in_mat}/SampleNumber_{N}.RData'
    output:
        heatmap = 'output/WTS4/normalize_1/{range}/{dist}/{in_mat}/SampleNumber_{N}/{v_type}/{l_comb}.png'
    params:
        label = 'data/WTS4_Eval_behavior_ACF.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{in_mat}/SampleNumber_{N}/{v_type}/{l_comb}.txt'
    conda:
        '../envs/myenv_WTS4_{dist}.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{in_mat}/SampleNumber_{N}/{v_type}/{l_comb}.log'
    shell:
        'src/WTS4_yshift_visualize.sh {input.RData} {params.label} {wildcards.v_type} {wildcards.l_comb} {output.heatmap} >& {log}'
