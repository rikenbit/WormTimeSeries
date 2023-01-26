# WTS5_SBD_abs
###################################################
N_SAMPLES = ["1"]
dist_data = ["SBD_abs"]
time_range = ["stimAfter"]
CellType = ["RIMR"]

rule all:
    input:
        # expand('output/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.RData',
        #     N=N_SAMPLES,
        #     dist=dist_data,
        #     celltype=CellType,
        #     range=time_range
        #     )
        expand('output/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/value_{celltype}.RData',
            N=N_SAMPLES,
            dist=dist_data,
            celltype=CellType,
            range=time_range
            )

rule WTS5_SBD_abs:
    input:
        'data/normalize_1/ReadData_{N}.RData'
    output:
        'output/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.RData'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.txt'
    conda:
        '../envs/myenv_WTS3_SBD.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.log'
    shell:
        'src/WTS5_SBD_abs.sh {wildcards.N} {input} {wildcards.range} {params.stim_xlsx} {wildcards.celltype} {output}  >& {log}'

rule yshift:
    input:
        'output/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_{celltype}.RData'
    output:
        'output/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/value_{celltype}.RData',
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/value_{celltype}.txt'
    conda:
        '../envs/myenv_WTS3_yshift.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS5/normalize_1/{range}/{dist}/SampleNumber_{N}/value_{celltype}.log'
    shell:
        'src/WTS5_yshift.sh {wildcards.N} {input} {params.stim_xlsx} {output} >& {log}'