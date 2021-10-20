# WTS4 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))

# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')


# Distance Data
dist_data = ["EUCL","SBD_abs"]
# data time range
time_range = ["all","stimAfter"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{dist}.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            )
        
rule DistMatrix:
    input:
        RData = 'data/normalize_1/ReadData_{N}.RData'
    output:
        dist_matrix = 'output/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{dist}.RData'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{dist}.txt'
    conda:
        '../envs/myenv_WTS4_{dist}.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{dist}.log'
    shell:
        'src/WTS4_{wildcards.dist}.sh {wildcards.N} {input.RData} {wildcards.range} {params.stim_xlsx} {output.dist_matrix}>& {log}'