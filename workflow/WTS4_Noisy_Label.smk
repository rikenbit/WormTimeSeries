# WTS4_Noisy_Label
###################################################
# N_SAMPLES = list(map(str, range(1, 29)))
# # remove artifact
# N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')
N_SAMPLES = ["2"]

NOISE_TEST = ["n1_noise"]

N_TRY = list(map(str, range(1, 11)))
# N_TRY = ["1"]

rule all:
    input:
        expand('data/{NOISE}_{TRY}/ReadData_{N}.RData',
            N=N_SAMPLES,
            NOISE=NOISE_TEST,
            TRY=N_TRY
            )
        
rule WTS4_Noisy_Label:
    input:
        ReadData = 'data/normalize_1/ReadData_{N}.RData'
    output:
        'data/{NOISE}_{TRY}/ReadData_{N}.RData'
    params:
        merged_cls = 'output/WTS4/normalize_1/stimAfter/SBD_abs/MCMIHOOI/Merged_cls/k_Number_9.RData'
    benchmark:
        'benchmarks/WTS4/{NOISE}_{TRY}/ReadData_{N}.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
        # tidyverseが入っているコンテナとして選んだ
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{NOISE}_{TRY}/ReadData_{N}.log'
    shell:
        'src/WTS4_Noisy_Label.sh {input.ReadData} {wildcards.N} {params.merged_cls} {wildcards.TRY} {output}>& {log}'