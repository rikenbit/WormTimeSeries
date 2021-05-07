# WTS3 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
### test ####
N_SAMPLES = N_SAMPLES[0]
# N_SAMPLES = N_SAMPLES[:4]
########
rule all:
    input:
        expand('output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.png', N=N_SAMPLES),
        expand('output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.RData', N=N_SAMPLES)

rule SBD:
    output:
        png = 'output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.png',
        RData = 'output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.RData'
    benchmark:
        'benchmarks/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.log'
    shell:
        'src/WTS3_SBD.sh {wildcards.N} {output.png} {output.RData} >& {log}'
###################################################