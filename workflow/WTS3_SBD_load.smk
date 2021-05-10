# WTS3 SBD load
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

rule all:
    input:
        expand('output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD_load.png', N=N_SAMPLES)

rule SBD_load:
    input:
        RData = 'output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.RData'
    output:
        png = 'output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD_load.png'
    benchmark:
        'benchmarks/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD_load.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD_load.log'
    shell:
        'src/WTS3_SBD_load.sh {wildcards.N} {output.png} {input.RData} >& {log}'
###################################################