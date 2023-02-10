# WTS5_tsPlot
###################################################
# N_SAMPLES = ["2"]
N_SAMPLES = ["1"]

CellType = [""]
# Shift_CellType = [""]
Shift_CellType = ["",""]

rule all:
    input:
        expand('output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_{N}_{celltype}_{s_celltype}.eps',
            N=N_SAMPLES,
            celltype=CellType,
            s_celltype=Shift_CellType
            )

rule WTS5_tsPlot:
    input:
        RData = 'data/normalize_1/ReadData_{N}.RData',
        stim = 'data/stimulation/stim_{N}.RData',
        # yshift ='output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_{N}/yshift_{celltype}.RData',
        # yshift_value = 'output/WTS3/normalize_1/stimAfter/SBD_abs_manual/SampleNumber_{N}/yshift_value_{celltype}.RData'
        yshift ='output/WTS5/normalize_1/stimAfter/SBD_abs/SampleNumber_{N}/yshift_{celltype}.RData',
        value_shift = 'output/WTS5/normalize_1/stimAfter/SBD_abs/SampleNumber_{N}/value_{celltype}.RData'
    output:
        'output/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_{N}_{celltype}_{s_celltype}.eps'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_{N}_{celltype}_{s_celltype}.txt'
    container:
        "docker://yamaken37/ggplot_svg:20230118"
    resources:
        mem_gb=200
    log:
        'logs/WTS5/normalize_1/stimAfter/SBD_abs/tsPlot/SampleNumber_{N}_{celltype}_{s_celltype}.log'
    shell:
        'src/WTS5_tsPlot.sh {wildcards.N} {wildcards.celltype} {wildcards.s_celltype} {input} {params.stim_xlsx} {output} >& {log}'