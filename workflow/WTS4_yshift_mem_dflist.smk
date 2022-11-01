# WTS4_yshift_mem_dflist
###################################################
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["6"]

rule all:
    input:
        expand('output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/DFs.RData',
            N_cls=N_CLUSTERS
            )
        
rule WTS4_yshift_mem_dflist:
    output:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/DFs.RData'
    params:
        sample_path = 'output/WTS4/normalize_1/stimAfter/SBD_abs/Distance'
    benchmark:
        'benchmarks/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/DFs.txt'
    container:
        "docker://yamaken37/dist_filter:20220624"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/DFs.log'
    shell:
        'src/WTS4_yshift_mem_dflist.sh {params.sample_path} {output} {wildcards.N_cls} >& {log}'        