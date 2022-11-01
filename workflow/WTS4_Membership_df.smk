# WTS4_Membership_df
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["6"]

rule all:
    input:
        expand('output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/SampleNumber_{N}.RData',
            N=N_SAMPLES,
            N_cls=N_CLUSTERS
            )
        
rule WTS4_Membership_df:
    input:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_F/k_Number_{N_cls}.RData'
    output:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/SampleNumber_{N}.RData'
    params:
        sample_path = 'output/WTS4/normalize_1/stimAfter/SBD_abs/Distance'
    benchmark:
        'benchmarks/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/SampleNumber_{N}.txt'
    container:
        "docker://yamaken37/dist_filter:20220624"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/SampleNumber_{N}.log'
    shell:
        'src/WTS4_Membership_df.sh {input} {output} {wildcards.N} {params.sample_path} >& {log}'        