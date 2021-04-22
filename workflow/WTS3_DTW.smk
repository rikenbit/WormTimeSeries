# WTS3 DTW 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

rule all:
    input:
        expand('output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW.png', N=N_SAMPLES),
        expand('output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW.RData', N=N_SAMPLES)

rule DTW:
    output:
        png = 'output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW.png',
        RData = 'output/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW.RData'
    benchmark:
        'benchmarks/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/DTW/normalize_1/all/SampleNumber_{N}/DTW.log'
    shell:
        'src/WTS3_DTW.sh {wildcards.N} {output.png} {output.RData} >& {log}'
###################################################