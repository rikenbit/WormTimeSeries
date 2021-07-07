# WTS3 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
# remove not annotation ASER, not periodic ASER →16sample
N_SAMPLES.remove('4')
N_SAMPLES.remove('7')
N_SAMPLES.remove('9')
N_SAMPLES.remove('15')
N_SAMPLES.remove('16')
N_SAMPLES.remove('22')
N_SAMPLES.remove('26')
N_SAMPLES.remove('28')

# Distance Data
dist_data = ["SBD"]
# data time range
time_range = ["all"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/SBD.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            ),
        expand('output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            )
        
rule SBD:
    input:
        RData = 'data/normalize_1/ReadData_{N}.RData'
    output:
        SBD = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/SBD.RData',
        yshift = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift.RData'
    params:
        args_shift = 'ASER'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/SBD.txt'
    conda:
        '../envs/myenv_WTS3_SBD.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/SBD.log'
    shell:
        'src/WTS3_SBD.sh {wildcards.N} {input.RData} {output.SBD} {output.yshift} {params.args_shift} >& {log}'
###################################################