# WTS4_DTW 
###################################################
# N_SAMPLES = list(map(str, range(1, 29)))

# # remove artifact
# N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')
N_SAMPLES = ["1"]


# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["DTW"]
# data time range
# time_range = ["all","stimAfter"]
time_range = ["stimAfter"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            )
        
rule WTS4_DTW:
    input:
        RData = 'data/normalize_1/ReadData_{N}.RData'
    output:
        dist_matrix = 'output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.txt'
    container:
        "docker://yamaken37/dtw:20220930"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.log'
    shell:
        'src/WTS4_DTW.sh {wildcards.N} {input.RData} {wildcards.range} {params.stim_xlsx} {output.dist_matrix}>& {log}'        