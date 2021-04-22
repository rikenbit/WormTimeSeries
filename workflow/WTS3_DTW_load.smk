# WTS3 DTW 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
### test ####
N_SAMPLES = N_SAMPLES[:4]
#######
rule all:
    input:
        expand('output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW_load.png', N=N_SAMPLES)

rule DTW:
    input:
        RData = 'output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW.RData'
    output:
        png = 'output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW_load.png'
    benchmark:
        'benchmarks/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW_load.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW_load.log'
    shell:
        'src/WTS3_DTW_load.sh {wildcards.N} {output.png} {input.RData} >& {log}'
###################################################