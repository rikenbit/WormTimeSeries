# WTS3 SBD load
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# option perplexity
op_per = ["15"]
# option max_iter
# op_max_iter = ["3000", "5000"]

rule all:
    input:
        expand('output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD_Purity_{op1}.png', N=N_SAMPLES, op1=op_per)

rule SBD_Purity:
    input:
        RData = 'output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.RData'
    output:
        png = 'output/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD_Purity_{op1}.png'
    benchmark:
        'benchmarks/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD_Purity_{op1}.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD_Purity_{op1}.log'
    shell:
        'src/WTS3_SBD_Purity.sh {wildcards.N} {output.png} {input.RData} {wildcards.op1} >& {log}'
###################################################