# WTS3 SBD
###################################################
N_SAMPLES = ["1"]
dist_data = ["SBD_abs"]
time_range = ["stimAfter"]
CellType = ["RIMR"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.RData',
            N=N_SAMPLES,
            dist=dist_data,
            celltype=CellType,
            range=time_range
            )

rule WTS5_SBD_abs:
    input:
        'data/normalize_1/ReadData_{N}.RData'
    output:
        yshift = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.RData'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx',
        args_shift = 'RIMR'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.txt'
    conda:
        '../envs/myenv_WTS3_SBD.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.log'
    shell:
        'src/WTS5_SBD_abs.sh {wildcards.N} {input} {wildcards.range} {params.stim_xlsx} {wildcards.celltype} {output}  >& {log}'
###################################################