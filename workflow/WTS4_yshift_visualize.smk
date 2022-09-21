# WTS4_yshift_visualize
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
# N_SAMPLES = ["21"]

# Distance Data
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# input matrix
# input_mat = ["Shift_F","Shift_F_M_ave","Shift_F_M_sd"]
input_matrix = ["Shift_F"]

# value type
# value_type =["zahlen","abs"]
value_type =["zahlen"]

# filter (label combination)
label_comb =["No_F","ALL","NaCl","1n","1p","1np"]
#### filter description
# "No_F" is no filter
# "ALL" = c("NaCl", "PC1_neg", "PC1_pos", "PC2", "PC3"),
# "NaCl" = c("NaCl"),
# "1n" = c("PC1_neg"),
# "1p" = c("PC1_pos"),
# "1np" = c("PC1_neg", "PC1_pos")
####

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
    container:
        "docker://yamaken37/yshift_visualize:20220921"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{in_mat}/SampleNumber_{N}/{v_type}/{l_comb}.log'
    shell:
        'src/WTS4_yshift_visualize.sh {input.RData} {params.label} {wildcards.v_type} {wildcards.l_comb} {output.heatmap} >& {log}'
