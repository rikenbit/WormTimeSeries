# WTS4_yshift_df
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

rule all:
    input:
        expand('output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_df/SampleNumber_{N}.RData',
            N=N_SAMPLES
            )
        
rule WTS4_yshift_df:
    input:
        yshift = 'output/WTS4/normalize_1/stimAfter/SBD_abs/Shift_F/SampleNumber_{N}.RData',
        dist = 'output/WTS4/normalize_1/stimAfter/SBD_abs/Distance/SampleNumber_{N}.RData'
    output:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_df/SampleNumber_{N}.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/stimAfter/SBD_abs/yShift_df/SampleNumber_{N}.txt'
    container:
        "docker://yamaken37/dist_filter:20220624"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/stimAfter/SBD_abs/yShift_df/SampleNumber_{N}.log'
    shell:
        'src/WTS4_yshift_df.sh {input.yshift} {input.dist} {output} >& {log}'