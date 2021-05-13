# WTS3 DTW 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# option perplexity
op_per = ["15"]

rule all:
    input:
        expand('output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW_load_{op1}.png', N=N_SAMPLES, op1=op_per)

rule SBD_load:
    input:
        RData = 'output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW.RData'
    output:
        png = 'output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW_load_{op1}.png'
    benchmark:
        'benchmarks/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW_load_{op1}.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW_load_{op1}.log'
    shell:
        'src/WTS3_DTW_load.sh {wildcards.N} {output.png} {input.RData} {wildcards.op1} >& {log}'
###################################################