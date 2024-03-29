# WTS3 DTW
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
dist_data = ["DTW"]
# data time range
time_range = ["all","stimAfter"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/DTW.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            )
        
rule DTW:
    input:
        RData = 'data/normalize_1/ReadData_{N}.RData'
    output:
        DTW = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/DTW.RData'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/DTW.txt'
    conda:
        '../envs/myenv_WTS3_DTW.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/DTW.log'
    shell:
        'src/WTS3_DTW.sh {wildcards.N} {input.RData} {wildcards.range} {params.stim_xlsx} {output.DTW}>& {log}'
###################################################