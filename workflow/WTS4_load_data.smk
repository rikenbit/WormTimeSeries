# WTS4_load_data
###################################################
N_SAMPLES = ["3","8","20","25"]
NOISE_TEST = ["n1_28sample"]

rule all:
    input:
        expand('data/{NOISE}/ReadData_{N}.RData',
            N=N_SAMPLES,
            NOISE=NOISE_TEST
            )
        
rule WTS4_load_data:
    output:
        'data/{NOISE}/ReadData_{N}.RData'
    params:
        animalname = 'data/normalize_1/AnimalName.csv'
    benchmark:
        'benchmarks/WTS4/{NOISE}/ReadData_{N}.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
        # tidyverseが入っているコンテナとして選んだ
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{NOISE}/ReadData_{N}.log'
    shell:
        'src/WTS4_load_data.sh {wildcards.N} {params.animalname} {output}>& {log}'