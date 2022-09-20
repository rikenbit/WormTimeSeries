# WTS4_yshift_visualize
###################################################
N_SAMPLES = list(map(str, range(1, 29)))

# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')


# Distance Data
dist_data = ["SBD_abs"]
# data time range
time_range = ["stimAfter"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Shift/SampleNumber_{N}.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/Shift_F/SampleNumber_{N}.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range
            )
        
rule WTS4_yshift_visualize:
    input:
        RData = 'data/normalize_1/ReadData_{N}.RData'
    output:
        Shift = 'output/WTS4/normalize_1/{range}/{dist}/Shift/SampleNumber_{N}.RData',
        Shift_F = 'output/WTS4/normalize_1/{range}/{dist}/Shift_F/SampleNumber_{N}.RData'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Shift/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS4_{dist}.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Shift/SampleNumber_{N}.log'
    shell:
        'src/WTS4_yshift_visualize.sh {wildcards.N} {input.RData} {wildcards.range} {params.stim_xlsx} {output.Shift} {output.Shift_F} >& {log}'
