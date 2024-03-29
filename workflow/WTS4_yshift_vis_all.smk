# WTS4_yshift_vis_all
###################################################
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["6"]

rule all:
    input:
        expand('output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis_all/k_Number_{N_cls}.png',
            N_cls=N_CLUSTERS
            )
        
rule WTS4_yshift_vis_all:
    input:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/Membership_df/k_Number_{N_cls}/DFs.RData'
    output:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis_all/k_Number_{N_cls}.png'
    params:
        sample_path = 'output/WTS4/normalize_1/stimAfter/SBD_abs/Distance'
    benchmark:
        'benchmarks/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis_all/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/dist_filter:20220624"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis_all/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_yshift_vis_all.sh {input} {output} >& {log}'        