# # WTS1 load raw
# ###################################################
# N_SAMPLES = list(map(str, range(1, 29)))

# rule all:
#     input:
#         expand('data/raw_CFP/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/raw_CFP/WTS1_sample_sheet.csv'),
#         expand('data/raw_CFP/AnimalName.csv'),
#         expand('data/raw_YFP/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/raw_YFP/WTS1_sample_sheet.csv'),
#         expand('data/raw_YFP/AnimalName.csv'),
#         expand('data/mCherry/mCherry_{N}.RData', N=N_SAMPLES),
#         expand('data/Position/Position_{N}.RData', N=N_SAMPLES),
#         expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)

# rule WTS1_load_raw_CFP:
#     output:
#         expand('data/raw_CFP/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/raw_CFP/WTS1_sample_sheet.csv'),
#         expand('data/raw_CFP/AnimalName.csv'),
#     params:
#         args1 = "raw_CFP",
#         args2 = "pi_k_Ch2"
#     benchmark:
#         'benchmarks/WTS1/load/raw_CFP.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/load/raw_CFP.log'
#     shell:
#         'src/WTS1_load_raw_n.sh {params.args1} {params.args2} >& {log}'

# # WTS1 load raw YFP
# rule WTS1_load_raw_YFP:
#     output:
#         expand('data/raw_YFP/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/raw_YFP/WTS1_sample_sheet.csv'),
#         expand('data/raw_YFP/AnimalName.csv'),
#     params:
#         args1 = "raw_YFP",
#         args2 = "pi_k_Ch3"
#     benchmark:
#         'benchmarks/WTS1/load/raw_YFP.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/load/raw_YFP.log'
#     shell:
#         'src/WTS1_load_raw_n.sh {params.args1} {params.args2} >& {log}'

# # WTS1 load raw other
# rule WTS1_load_raw_other:
#     output:
#         expand('data/mCherry/mCherry_{N}.RData', N=N_SAMPLES),
#         expand('data/Position/Position_{N}.RData', N=N_SAMPLES),
#         expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)
#     benchmark:
#         'benchmarks/WTS1/load/raw_other.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/load/raw_other.log'
#     shell:
#         'src/WTS1_load_raw_other.sh >& {log}'
# ###################################################

# # WTS1 load normalize
# ###################################################
# N_SAMPLES = list(map(str, range(1, 29)))
# rule all:
#     input:
#         expand('data/normalize_1/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/normalize_1/WTS1_sample_sheet.csv'),
#         expand('data/normalize_1/AnimalName.csv'),
#         expand('data/normalize_2/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/normalize_2/WTS1_sample_sheet.csv'),
#         expand('data/normalize_2/AnimalName.csv'),
#         expand('data/normalize_3/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/normalize_3/WTS1_sample_sheet.csv'),
#         expand('data/normalize_3/AnimalName.csv'),
#         expand('data/normalize_4/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/normalize_4/WTS1_sample_sheet.csv'),
#         expand('data/normalize_4/AnimalName.csv')

# rule WTS1_load_normalize_1:
#     output:
#         expand('data/normalize_1/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/normalize_1/WTS1_sample_sheet.csv'),
#         expand('data/normalize_1/AnimalName.csv'),
#     params:
#         args1 = "normalize_1"
#     benchmark:
#         'benchmarks/WTS1/load/normalize_1.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/load/normalize_1.log'
#     shell:
#         'src/WTS1_load_normalize_n.sh {params.args1} >& {log}'
# rule WTS1_load_normalize_2:
#     output:
#         expand('data/normalize_2/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/normalize_2/WTS1_sample_sheet.csv'),
#         expand('data/normalize_2/AnimalName.csv'),
#     params:
#         args1 = "normalize_2"
#     benchmark:
#         'benchmarks/WTS1/load/normalize_2.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/load/normalize_2.log'
#     shell:
#         'src/WTS1_load_normalize_n.sh {params.args1} >& {log}'
# rule WTS1_load_normalize_3:
#     output:
#         expand('data/normalize_3/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/normalize_3/WTS1_sample_sheet.csv'),
#         expand('data/normalize_3/AnimalName.csv'),
#     params:
#         args1 = "normalize_3"
#     benchmark:
#         'benchmarks/WTS1/load/normalize_3.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/load/normalize_3.log'
#     shell:
#         'src/WTS1_load_normalize_n.sh {params.args1} >& {log}'
# rule WTS1_load_normalize_4:
#     output:
#         expand('data/normalize_4/ReadData_{N}.RData', N=N_SAMPLES),
#         expand('data/normalize_4/WTS1_sample_sheet.csv'),
#         expand('data/normalize_4/AnimalName.csv'),
#     params:
#         args1 = "normalize_4"
#     benchmark:
#         'benchmarks/WTS1/load/normalize_4.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/load/normalize_4.log'
#     shell:
#         'src/WTS1_load_normalize_n.sh {params.args1} >& {log}'
# ###################################################

# # WTS1 plot CFP
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"
# DATA_DIR = config['DATA_DIR_CFP']
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET_CFP'], dtype='string')

# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

# rule WTS1_plot_CFP:
#     input:
#         expand('data/{DATA_DIR}/ReadData_{Sample}.RData', Sample = paramspace["SampleNumber"], DATA_DIR = DATA_DIR),
#         expand('data/stimulation/stim_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/Position/Position_{Sample}.RData', Sample = paramspace["SampleNumber"])
#     output:
#         f"output/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = {DATA_DIR}
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_plot.sh {params.args1} {params.args2} {params.args3} {params.args4} {output} >& {log}'
# ###################################################

# # WTS1 plot YFP
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"
# DATA_DIR = config['DATA_DIR_YFP']
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET_YFP'], dtype='string')

# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

# rule WTS1_plot_YFP:
#     input:
#         expand('data/{DATA_DIR}/ReadData_{Sample}.RData', Sample = paramspace["SampleNumber"], DATA_DIR = DATA_DIR),
#         expand('data/stimulation/stim_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/Position/Position_{Sample}.RData', Sample = paramspace["SampleNumber"])
#     output:
#         f"output/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = {DATA_DIR}
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_plot.sh {params.args1} {params.args2} {params.args3} {params.args4} {output} >& {log}'
# ###################################################

# # WTS1 plot normalize_1
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"
# DATA_DIR = config['DATA_DIR_normalize_1']
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET_normalize_1'], dtype='string')

# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

# rule WTS1_plot_normalize_1:
#     input:
#         expand('data/{DATA_DIR}/ReadData_{Sample}.RData', Sample = paramspace["SampleNumber"], DATA_DIR = DATA_DIR),
#         expand('data/stimulation/stim_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/Position/Position_{Sample}.RData', Sample = paramspace["SampleNumber"])
#     output:
#         f"output/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = {DATA_DIR}
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_plot.sh {params.args1} {params.args2} {params.args3} {params.args4} {output} >& {log}'
# ###################################################

# # WTS1 plot normalize_2
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"
# DATA_DIR = config['DATA_DIR_normalize_2']
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET_normalize_2'], dtype='string')

# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

# rule WTS1_plot_normalize_2:
#     input:
#         expand('data/{DATA_DIR}/ReadData_{Sample}.RData', Sample = paramspace["SampleNumber"], DATA_DIR = DATA_DIR),
#         expand('data/stimulation/stim_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/Position/Position_{Sample}.RData', Sample = paramspace["SampleNumber"])
#     output:
#         f"output/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = {DATA_DIR}
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_plot.sh {params.args1} {params.args2} {params.args3} {params.args4} {output} >& {log}'
# ###################################################

# # WTS1 plot normalize_3
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"
# DATA_DIR = config['DATA_DIR_normalize_3']
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET_normalize_3'], dtype='string')

# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

# rule WTS1_plot_normalize_3:
#     input:
#         expand('data/{DATA_DIR}/ReadData_{Sample}.RData', Sample = paramspace["SampleNumber"], DATA_DIR = DATA_DIR),
#         expand('data/stimulation/stim_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/Position/Position_{Sample}.RData', Sample = paramspace["SampleNumber"])
#     output:
#         f"output/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = {DATA_DIR}
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_plot.sh {params.args1} {params.args2} {params.args3} {params.args4} {output} >& {log}'
# ###################################################

# # WTS1 plot normalize_4
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"
# DATA_DIR = config['DATA_DIR_normalize_4']
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET_normalize_4'], dtype='string')

# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

# rule WTS1_plot_normalize_4:
#     input:
#         expand('data/{DATA_DIR}/ReadData_{Sample}.RData', Sample = paramspace["SampleNumber"], DATA_DIR = DATA_DIR),
#         expand('data/stimulation/stim_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/mCherry/mCherry_{Sample}.RData', Sample = paramspace["SampleNumber"]),
#         expand('data/Position/Position_{Sample}.RData', Sample = paramspace["SampleNumber"])
#     output:
#         f"output/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = {DATA_DIR}
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_plot.sh {params.args1} {params.args2} {params.args3} {params.args4} {output} >& {log}'
# ###################################################

# # WTS1 heatmap
# ###################################################
# DATA_DIR = ["raw_CFP", "raw_YFP", "normalize_1", "normalize_2", "normalize_3", "normalize_4"]
# SAMPLES = list(map(str, range(1, 29)))

# rule all:
#     input:
#         expand('output/WTS1/heatmap/{DATA_DIR}/SampleNumber_{SAMPLES}.png', DATA_DIR = DATA_DIR, SAMPLES = SAMPLES)

# rule heatmap:
#     output:
#             'output/WTS1/heatmap/{DATA_DIR}/SampleNumber_{SAMPLES}.png'
#     benchmark:
#             'benchmarks/WTS1/heatmap/{DATA_DIR}/SampleNumber_{SAMPLES}.txt'
#     conda:
#             'envs/myenv_WTS1.yaml'
#     resources:
#             mem_gb=200
#     log:
#             'logs/WTS1/heatmap/{DATA_DIR}/SampleNumber_{SAMPLES}.log'
#     shell:
#             'src/WTS1_heatmap.sh {wildcards.SAMPLES} {wildcards.DATA_DIR} {output} >& {log}'
# ###################################################

# # WTS2 correlogram
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"
# # DATA = ["raw_CFP", "raw_YFP", "normalize_1", "normalize_2", "normalize_3", "normalize_4"]
# raw_CFP = pd.read_csv(config['raw_CFP'], dtype='string')
# raw_CFP['Data'] = 'raw_CFP'
# raw_YFP = pd.read_csv(config['raw_YFP'], dtype='string')
# raw_YFP['Data'] = 'raw_YFP'
# normalize_1 = pd.read_csv(config['normalize_1'], dtype='string')
# normalize_1['Data'] = 'normalize_1'
# normalize_2 = pd.read_csv(config['normalize_2'], dtype='string')
# normalize_2['Data'] = 'normalize_2'
# normalize_3 = pd.read_csv(config['normalize_3'], dtype='string')
# normalize_3['Data'] = 'normalize_3'
# normalize_4 = pd.read_csv(config['normalize_4'], dtype='string')
# normalize_4['Data'] = 'normalize_4'
# df_Data = pd.concat([raw_CFP, raw_YFP, normalize_1, normalize_2, normalize_3, normalize_4], axis=0)

# # TimeFrame = ["all", "before", "after"]
# TF_all = df_Data.copy()
# TF_all['TF'] = 'all'
# TF_before = df_Data.copy()
# TF_before['TF'] = 'before'
# TF_after = df_Data.copy()
# TF_after['TF'] = 'after'

# df_TF = pd.concat([TF_all,TF_before,TF_after], axis=0)

# # LAG_MAX = ["50", "100"]
# LAG_50 = df_TF.copy()
# LAG_50['LAG'] = '50'
# LAG_300 = df_TF.copy()
# LAG_300['LAG'] = '300'

# df_LAG = pd.concat([LAG_50,LAG_300], axis=0)
# # ACF = ["Acf","pAcf"]
# ACF_Acf = df_LAG.copy()
# ACF_Acf['ACF'] = 'Acf'
# ACF_pAcf = df_LAG.copy()
# ACF_pAcf['ACF'] = 'pAcf'

# df_ACF = pd.concat([ACF_Acf,ACF_pAcf], axis=0)
# df = df_ACF.reset_index(drop=True)
# df = df.reindex(columns=['Data', 'TF', 'LAG', 'ACF', 'SampleNumber','CellNumber', 'CellType'])
# # paramspace = Paramspace(df, filename_params=['CellNumber', 'CellType'], param_sep="_")
# # #### test####
# # df_test = df.copy()
# # df_test = df_test[df_test['Data'].isin(['normalize_1']) & df_test['TF'].isin(['before']) & df_test['LAG'].isin(['300']) & df_test['ACF'].isin(['Acf']) & df_test['SampleNumber'].isin(['13']) & df_test['CellType'].isin(['ASER'])]
# # paramspace = Paramspace(df_test, filename_params=['CellNumber', 'CellType'], param_sep="_")
# # #### test####
# ### filter####
# df_test = df.copy()
# df_test = df_test[df_test['Data'].isin(['normalize_1'])  & df_test['LAG'].isin(['300']) & ~df_test['SampleNumber'].isin(['3','8','20','25'])]
# paramspace = Paramspace(df_test, filename_params=['CellNumber', 'CellType'], param_sep="_")
# ### filter####

# rule all:
#     input:
#     	# output/WTS2/correlogram/raw_CFP/all/50/Acf/SampleNumber_1/CellNumber_87_CellType_ADAR.png
#         expand('output/WTS2/correlogram/{params}.png', params = paramspace.instance_patterns)
# rule correlogram:
#     output:
#     	# f"output/WTS2/correlogram/{DATA_DIR}/{TimeFrame}/{LAG_MAX}/{ACF}/{paramspace.wildcard_pattern}.png"
#     	expand('output/WTS2/correlogram/{params}.png', params = paramspace.wildcard_pattern)
#     params:
#     	args1 = lambda w: w["Data"],
#     	args2 = lambda w: w["TF"],
#     	args3 = lambda w: w["LAG"],
#     	args4 = lambda w: w["ACF"],
#         args5 = lambda w: w["SampleNumber"],
#         args6 = lambda w: w["CellNumber"],
#         args7 = lambda w: w["CellType"]
#     benchmark:
#         f'benchmarks/WTS2/correlogram/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS2.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS2/correlogram/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS2_correlogram.sh {params.args1} {params.args2} {params.args3} {params.args4} {params.args5} {params.args6} {params.args7} {output} >& {log}'
# ###################################################

# # WTS2 heatmap normalize_1
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"

# # DATA = ["raw_CFP", "raw_YFP", "normalize_1", "normalize_2", "normalize_3", "normalize_4"]
# # normalize_1 = pd.read_csv(config['normalize_1'], dtype='string')
# # normalize_1['Data'] = 'normalize_1'
# # normalize_2 = pd.read_csv(config['normalize_2'], dtype='string')
# # normalize_2['Data'] = 'normalize_2'
# # normalize_3 = pd.read_csv(config['normalize_3'], dtype='string')
# # normalize_3['Data'] = 'normalize_3'
# # normalize_4 = pd.read_csv(config['normalize_4'], dtype='string')
# # normalize_4['Data'] = 'normalize_4'
# # df_Data = pd.concat([raw_CFP, raw_YFP, normalize_1, normalize_2, normalize_3, normalize_4], axis=0)
# df_Data = pd.read_csv(config['normalize_1'], dtype='string')
# df_Data['Data'] = 'normalize_1'

# # TimeFrame = ["all", "before", "after"]
# TF_all = df_Data.copy()
# TF_all['TF'] = 'all'
# TF_before = df_Data.copy()
# TF_before['TF'] = 'before'
# TF_after = df_Data.copy()
# TF_after['TF'] = 'after'

# df_TF = pd.concat([TF_all,TF_before,TF_after], axis=0)

# # LAG_MAX = ["50", "300"]
# LAG_1 = df_TF.copy()
# LAG_1['LAG'] = '1'
# df_LAG = LAG_1.copy()
# lag = list(map(str, range(2, 301)))
# for lags in lag:
# 	LAG_N = df_TF.copy()
# 	LAG_N['LAG'] = lags
# 	df_LAG = pd.concat([df_LAG,LAG_N], axis=0)
# #
# df = df_LAG.reset_index(drop=True)
# df = df.reindex(columns=['Data', 'TF', 'LAG', 'SampleNumber','CellNumber', 'CellType'])


# ### filter####
# df_test = df.copy()
# df_test = df_test[~df_test['SampleNumber'].isin(['3','8','20','25'])]

# #### test####
# df_test = df_test[df_test['TF'].isin(['all']) & df_test['LAG'].isin(['1']) & df_test['SampleNumber'].isin(['13'])]
# #### test####

# paramspace = Paramspace(df_test, filename_params=['CellNumber', 'CellType'], param_sep="_")
# ### filter####

# rule all:
#     input:
#     	# output/WTS2/correlogram/raw_CFP/all/50/Acf/SampleNumber_1/CellNumber_87_CellType_ADAR.png
#         expand('output/WTS2/correlogram/{params}.png', params = paramspace.instance_patterns)
# rule WTS2 heatmap:
#     output:
#     	# f"output/WTS2/correlogram/{DATA_DIR}/{TimeFrame}/{LAG_MAX}/{ACF}/{paramspace.wildcard_pattern}.png"
#     	expand('output/WTS2/correlogram/{params}.png', params = paramspace.wildcard_pattern)
#     params:
#     	args1 = lambda w: w["Data"],
#     	args2 = lambda w: w["TF"],
#     	args3 = lambda w: w["LAG"],
#     	args4 = lambda w: w["ACF"],
#         args5 = lambda w: w["SampleNumber"],
#         args6 = lambda w: w["CellNumber"],
#         args7 = lambda w: w["CellType"]
#     benchmark:
#         f'benchmarks/WTS2/correlogram/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS2.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS2/correlogram/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS2_correlogram.sh {params.args1} {params.args2} {params.args3} {params.args4} {params.args5} {params.args6} {params.args7} {output} >& {log}'
# ###################################################

# WTS2 heatmap normalize_1
###################################################
Data = 'normalize_1'
TF = ["all", "before", "after"]
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
LAG = list(map(str, range(1, 301)))
### filter####
TF = TF[0]
N_SAMPLES = N_SAMPLES[0]
LAG = LAG[0]
### filter####

rule all:
    input:
    	# output/WTS2/heatmap/normalize_1/all/SampleNumber_1/τ1.png
    	expand('output/WTS2/heatmap/{D}/{T}/SampleNumber_{N}/τ{L}.png',
    	 D=Data,
    	 T=TF,
    	 N=N_SAMPLES,
    	 L=LAG)
    	
rule heatmap:
    output:
    	'output/WTS2/heatmap/{D}/{T}/SampleNumber_{N}/τ{L}.png'
    benchmark:
        'benchmarks/WTS2/heatmap/{D}/{T}/SampleNumber_{N}/τ{L}.txt'
    conda:
        'envs/myenv_WTS2.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS2/heatmap/{D}/{T}/SampleNumber_{N}/τ{L}.log'
    shell:
        'src/WTS2_heatmap.sh {wildcards.D} {wildcards.T} {wildcards.L} {wildcards.N} {output} >& {log}'
###################################################