# WTS4_Membership_vis
###################################################
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["6"]

N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

rule all:
    input:
        expand('output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis/k_Number_{N_cls}/SampleNumber_{N}.png',
            N_cls=N_CLUSTERS,
            N=N_SAMPLES
            )
        
rule WTS4_Membership_vis:
    input:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/SampleNumber_{N}.RData'
    output:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis/k_Number_{N_cls}/SampleNumber_{N}.png'
    benchmark:
        'benchmarks/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis/k_Number_{N_cls}/SampleNumber_{N}.txt'
    container:
        "docker://yamaken37/dist_filter:20220624"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/stimAfter/SBD_abs/Membership_vis/k_Number_{N_cls}/SampleNumber_{N}.log'
    shell:
        'src/WTS4_Membership_vis.sh {input} {output} >& {log}'   