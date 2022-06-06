# WTS2_correlogram_cells
###################################################
import pandas as pd
from snakemake.utils import Paramspace

configfile: "config.yaml"
normalize_1 = pd.read_csv(config['normalize_1'], dtype='string')
normalize_1['Data'] = 'normalize_1'
df_Data = pd.concat([normalize_1], axis=0)

# TimeFrame = ["all", "before", "after"]
TF_all = df_Data.copy()
TF_all['TF'] = 'all'
TF_before = df_Data.copy()
TF_before['TF'] = 'before'
TF_after = df_Data.copy()
TF_after['TF'] = 'after'

df_TF = pd.concat([TF_all,TF_before,TF_after], axis=0)

# LAG_MAX = ["50", "100"]
LAG_50 = df_TF.copy()
LAG_50['LAG'] = '50'
LAG_300 = df_TF.copy()
LAG_300['LAG'] = '600'

df_LAG = pd.concat([LAG_50,LAG_300], axis=0)
# ACF = ["Acf","pAcf"]
ACF_Acf = df_LAG.copy()
ACF_Acf['ACF'] = 'Acf'
ACF_pAcf = df_LAG.copy()
ACF_pAcf['ACF'] = 'pAcf'

df_ACF = pd.concat([ACF_Acf,ACF_pAcf], axis=0)
df = df_ACF.reset_index(drop=True)
df = df.reindex(columns=['TF', 'LAG', 'ACF', 'SampleNumber', 'CellType'])
### filter####
df_test = df.copy()
# df_test = df_test[df_test['Data'].isin(['normalize_1']) & df_test['LAG'].isin(['600']) & ~df_test['SampleNumber'].isin(['3','8','20','25']) & df_test['CellType'].isin(['ASER']) & df_test['ACF'].isin(['Acf']) & df_test['TF'].isin(['after'])]
df_test = df_test[df_test['LAG'].isin(['600']) & df_test['ACF'].isin(['Acf']) & df_test['TF'].isin(['after'])]
paramspace = Paramspace(df_test, filename_params=['SampleNumber'], param_sep="_")
### filter####

rule all:
    input:
        # output/WTS2/correlogram/normalize_1/raw_CFP/all/50/Acf/SampleNumber_1/CellNumber_87_CellType_ADAR.png
        expand('output/WTS2/correlogram/normalize_1/{params}.png', params = paramspace.instance_patterns)
rule WTS2_correlogram_cells:
    output:
        # f"output/WTS2/correlogram/normalize_1/{DATA_DIR}/{TimeFrame}/{LAG_MAX}/{ACF}/{paramspace.wildcard_pattern}.png"
        expand('output/WTS2/correlogram/normalize_1/{params}.png', params = paramspace.wildcard_pattern)
    params:
        # args1 = lambda w: w["Data"],
        args2 = lambda w: w["TF"],
        args3 = lambda w: w["LAG"],
        args4 = lambda w: w["ACF"],
        args5 = lambda w: w["SampleNumber"],
        # args6 = lambda w: w["CellNumber"],
        args7 = lambda w: w["CellType"]
    benchmark:
        f'benchmarks/WTS2/correlogram/normalize_1/{paramspace.wildcard_pattern}.txt'
    conda:
        'envs/myenv_WTS2.yaml'
    resources:
        mem_gb=200
    log:
        f'logs/WTS2/correlogram/normalize_1/{paramspace.wildcard_pattern}.log'
    shell:
        'src/WTS2_correlogram.sh {params.args2} {params.args3} {params.args4} {params.args5} {params.args7} {output} >& {log}'