# WTS3 EUCL 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

rule all:
    input:
        expand('output/WTS3/EUCL/normalize_1/all/SampleNumber_{N}/EUCL.png', N=N_SAMPLES)

rule EUCL:
    output:
        'output/WTS3/EUCL/normalize_1/all/SampleNumber_{N}/EUCL.png'
    benchmark:
        'benchmarks/WTS3/EUCL/normalize_1/all/SampleNumber_{N}/EUCL.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/EUCL/normalize_1/all/SampleNumber_{N}/EUCL.log'
    shell:
        'src/WTS3_EUCL.sh {wildcards.N} {output} >& {log}'
###################################################