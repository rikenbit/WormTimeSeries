# WTS3 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
# remove not annotation ASER, not periodic ASER →17sample
N_SAMPLES.remove('4')
# N_SAMPLES.remove('7')
N_SAMPLES.remove('9')
N_SAMPLES.remove('15')
N_SAMPLES.remove('16')
N_SAMPLES.remove('22')
N_SAMPLES.remove('26')
N_SAMPLES.remove('28')

# Distance Data
dist_data = ["SBD","SBD_abs"]
# data time range
# time_range = ["all","stimAfter"]
time_range = ["stimAfter"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_value.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            )
rule yshift:
    input:
        RData = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift.RData'
    output:
        yshift_value = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_value.RData',
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_value.txt'
    conda:
        '../envs/myenv_WTS3_yshift.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_value.log'
    shell:
        'src/WTS3_yshift.sh {wildcards.N} {input.RData} {params.stim_xlsx} {output.yshift_value} >& {log}'
###################################################