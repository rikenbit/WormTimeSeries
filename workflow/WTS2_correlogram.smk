# WTS2 correlogram
###################################################
import pandas as pd
from snakemake.utils import Paramspace

configfile: "config.yaml"

raw_CFP = pd.read_csv(config['raw_CFP'], dtype='string')
raw_CFP['Data'] = 'raw_CFP'
raw_YFP = pd.read_csv(config['raw_YFP'], dtype='string')
raw_YFP['Data'] = 'raw_YFP'
normalize_1 = pd.read_csv(config['normalize_1'], dtype='string')
normalize_1['Data'] = 'normalize_1'
normalize_2 = pd.read_csv(config['normalize_2'], dtype='string')
normalize_2['Data'] = 'normalize_2'
normalize_3 = pd.read_csv(config['normalize_3'], dtype='string')
normalize_3['Data'] = 'normalize_3'
normalize_4 = pd.read_csv(config['normalize_4'], dtype='string')
normalize_4['Data'] = 'normalize_4'
df_Data = pd.concat([raw_CFP, raw_YFP, normalize_1, normalize_2, normalize_3, normalize_4], axis=0)

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
LAG_300['LAG'] = '300'
LAG_500 = df_TF.copy()
LAG_500['LAG'] = '600'

df_LAG = pd.concat([LAG_50,LAG_300,LAG_500], axis=0)
# ACF = ["Acf","pAcf"]
ACF_Acf = df_LAG.copy()
ACF_Acf['ACF'] = 'Acf'
ACF_pAcf = df_LAG.copy()
ACF_pAcf['ACF'] = 'pAcf'

df_ACF = pd.concat([ACF_Acf,ACF_pAcf], axis=0)
df = df_ACF.reset_index(drop=True)
df = df.reindex(columns=['Data', 'TF', 'LAG', 'ACF', 'SampleNumber','CellNumber', 'CellType'])
# paramspace = Paramspace(df, filename_params=['CellNumber', 'CellType'], param_sep="_")

### filter####
df_test = df.copy()
df_test = df_test[df_test['Data'].isin(['normalize_1']) & df_test['TF'].isin(['all']) & df_test['LAG'].isin(['500']) & df_test['ACF'].isin(['Acf']) & ~df_test['SampleNumber'].isin(['3','8','20','25'])]
paramspace = Paramspace(df_test, filename_params=['CellNumber', 'CellType'], param_sep="_")
### filter####

rule all:
    input:
        expand('output/WTS2/correlogram/{params}.png', params = paramspace.instance_patterns)

rule correlogram:
    output:
        expand('output/WTS2/correlogram/{params}.png', params = paramspace.wildcard_pattern)
    params:
        args1 = lambda w: w["Data"],
        args2 = lambda w: w["TF"],
        args3 = lambda w: w["LAG"],
        args4 = lambda w: w["ACF"],
        args5 = lambda w: w["SampleNumber"],
        args6 = lambda w: w["CellNumber"],
        args7 = lambda w: w["CellType"]
    benchmark:
        f'benchmarks/WTS2/correlogram/{paramspace.wildcard_pattern}.txt'
    conda:
        '../envs/myenv_WTS2.yaml'
    resources:
        mem_gb=200
    log:
        f'logs/WTS2/correlogram/{paramspace.wildcard_pattern}.log'
    shell:
        'src/WTS2_correlogram.sh {params.args1} {params.args2} {params.args3} {params.args4} {params.args5} {params.args6} {params.args7} {output} >& {log}'
###################################################