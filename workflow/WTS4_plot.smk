import pandas as pd
from snakemake.utils import Paramspace

DATA_DIR = ["n1_28sample"]
SAMPLE_SHEET = pd.read_csv('data/n1_28sample/WTS4_sample_sheet.csv', dtype='string')

paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

rule all:
    input:
        expand('output/WTS4/{DATA_DIR}/plot/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

rule WTS1_plot_normalize_1:
    input:
        expand('data/{DATA_DIR}/ReadData_{Sample}.RData', Sample = paramspace["SampleNumber"], DATA_DIR = DATA_DIR),
        expand('data/stimulation/stim_{Sample}.RData', Sample = paramspace["SampleNumber"]),
        expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
        expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
        expand('data/Position/Position_{Sample}.RData', Sample = paramspace["SampleNumber"])
    output:
        f"output/WTS4/{DATA_DIR}/plot/{paramspace.wildcard_pattern}.png"
    params:
        args1 = lambda w: w["SampleNumber"],
        args2 = lambda w: w["CellNumber"],
        args3 = lambda w: w["CellType"],
        args4 = {DATA_DIR}
    benchmark:
        f'benchmarks/WTS4/{DATA_DIR}/plot/{paramspace.wildcard_pattern}.txt'
    container:
        "docker://yamaken37/plot_tf:20220314"
    resources:
        mem_gb=200
    log:
        f'logs/WTS4/{DATA_DIR}/plot/{paramspace.wildcard_pattern}.log'
    shell:
        'src/WTS4_plot.sh {params.args1} {params.args2} {params.args3} {params.args4} {output} >& {log}'