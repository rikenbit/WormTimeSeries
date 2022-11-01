# WTS4_yshift_vis
###################################################

N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

rule all:
    input:
        expand('output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis/SampleNumber_{N}.png',
            N=N_SAMPLES
            )
        
rule WTS4_yshift_vis:
    input:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_df/SampleNumber_{N}.RData'
    output:
        'output/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis/SampleNumber_{N}.png'
    benchmark:
        'benchmarks/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis/SampleNumber_{N}.txt'
    container:
        "docker://yamaken37/dimreduc:20211213"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/stimAfter/SBD_abs/yShift_vis/SampleNumber_{N}.log'
    shell:
        'src/WTS4_yshift_vis.sh {input} {output} >& {log}'