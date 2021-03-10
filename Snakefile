# # WTS1 load raw CFP→エラー
# ###################################################
# # N_SAMPLES = list(map(str, range(1, 29)))
# N_SAMPLES = list(map(str, range(1, 2)))
# N_DATA_RAW = ["raw_CFP", "raw_YFP"]
# rule all:
# 	input:
# 		expand('data/{D}/ReadData_{N}.RData', N=N_SAMPLES, D=N_DATA_RAW),
#         expand('data/{D}/WTS1_sample_sheet.csv', D=N_DATA_RAW),
#         expand('data/{D}/AnimalName.csv', D=N_DATA_RAW),
#         expand('data/mCherry/mCherry_{N}.RData', N=N_SAMPLES),
#         expand('data/Position/Position_{N}.RData', N=N_SAMPLES),
#         expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)

# rule WTS1_load:
#     output:
#         expand('data/{D}/ReadData_{N}.RData', N=N_SAMPLES, D=N_DATA_RAW),
#         expand('data/{D}/WTS1_sample_sheet.csv', D=N_DATA_RAW),
#         expand('data/{D}/AnimalName.csv', D=N_DATA_RAW),
#         expand('data/mCherry/mCherry_{N}.RData', N=N_SAMPLES),
#         expand('data/Position/Position_{N}.RData', N=N_SAMPLES),
#         expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)
#     params:
#         args1 = N_DATA_RAW
#     benchmark:
#         'benchmarks/WTS1/{D}/WTS1_load.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/{D}/WTS1_load.log'
#     shell:
#         'src/WTS1_load.sh {wildcards.D} >& {log}'
# ###################################################

# WTS1 load raw CFP
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_DATA_RAW = ["raw_CFP"]
rule all:
    input:
        expand('data/raw_CFP/ReadData_{N}.RData', N=N_SAMPLES),
        expand('data/raw_CFP/WTS1_sample_sheet.csv'),
        expand('data/raw_CFP/AnimalName.csv'),
        expand('data/mCherry/mCherry_{N}.RData', N=N_SAMPLES),
        expand('data/Position/Position_{N}.RData', N=N_SAMPLES),
        expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)

rule WTS1_load:
    output:
        expand('data/raw_CFP/ReadData_{N}.RData', N=N_SAMPLES),
        expand('data/raw_CFP/WTS1_sample_sheet.csv'),
        expand('data/raw_CFP/AnimalName.csv'),
        expand('data/mCherry/mCherry_{N}.RData', N=N_SAMPLES),
        expand('data/Position/Position_{N}.RData', N=N_SAMPLES),
        expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)
    params:
        args1 = N_DATA_RAW
    benchmark:
        'benchmarks/WTS1/raw_CFP/WTS1_load.txt'
    conda:
        'envs/myenv_WTS1.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS1/raw_CFP/WTS1_load.log'
    shell:
        'src/WTS1_load.sh {params.args1} >& {log}'
###################################################

# # WTS1 load normalize
# ###################################################
# N_SAMPLES = list(map(str, range(1, 29)))
# N_DATA_N = ["normalize_1", "normalize_2", "normalize_3", "normalize_4"]
# rule all:
#     input:
#         expand('data/{D}/ReadData_{N}.RData', N=N_SAMPLES, D=N_DATA_N,
#         expand('data/{D}/WTS1_sample_sheet.csv', D=N_DATA_N,
#         expand('data/{D}/AnimalName.csv', D=N_DATA_N,
#         expand('data/mCherry/mCherry_{N}.RData', N=N_SAMPLES),
#         expand('data/Position/Position_{N}.RData', N=N_SAMPLES),
#         expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)

# rule WTS1_load:
#     output:
#         expand('data/{D}/ReadData_{N}.RData', N=N_SAMPLES, D=N_DATA_N,
#         expand('data/{D}/WTS1_sample_sheet.csv', D=N_DATA_N,
#         expand('data/{D}/AnimalName.csv', D=N_DATA_N,
#         expand('data/mCherry/mCherry_{N}.RData', N=N_SAMPLES),
#         expand('data/Position/Position_{N}.RData', N=N_SAMPLES),
#         expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)
#     benchmark:
#         'benchmarks/WTS1/{params.args1}/WTS1_load.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/{params.args1}/WTS1_load.log'
#     shell:
#         'src/WTS1_load.sh {params.args1} >& {log}'
# ###################################################

# # WTS1 plot raw_CFP
# ###################################################
# import pandas as pd
# # from snakemake.utils import min_version
# from snakemake.utils import Paramspace


# # min_version("5.3.2")
# configfile: "config/config_raw_CFP.yaml"

# # read sample_sheet
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')

# DATA_DIR = config['DATA_DIR']

# # paramspace
# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

# rule WTS1_1cell:
#     output:
#         f"output/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = DATA_DIR
#     benchmark:
#         f'benchmarks/WTS1/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}'
# ###################################################

# # WTS1 plot raw_YFP
# ###################################################
# import pandas as pd
# # from snakemake.utils import min_version
# from snakemake.utils import Paramspace


# # min_version("5.3.2")
# configfile: "config/config_raw_YFP.yaml"

# # read sample_sheet
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')

# DATA_DIR = config['DATA_DIR']

# # paramspace
# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace.instance_patterns, DATA_DIR = DATA_DIR)

# rule WTS1_1cell:
#     output:
#         f"output/WTS1/plot/{DATA_DIR}/{paramspace.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = DATA_DIR
#     benchmark:
#         f'benchmarks/WTS1/{paramspace.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/{paramspace.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}'
# ###################################################

# # WTS2 correlogram τ50
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"

# # read sample_sheet
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')
# # remove sample
# SAMPLE_SHEET = SAMPLE_SHEET[~SAMPLE_SHEET['SampleNumber'].isin(['7','12','13','14'])]
# SAMPLE_SHEET  = SAMPLE_SHEET .reset_index(drop=True)

# # paramspace
# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
# 	input:
# 		expand('output/WTS2/ACFτ50/{params}_Acf.png', params = paramspace.instance_patterns),
# 		expand('output/WTS2/ACFτ50/{params}_pAcf.png', params = paramspace.instance_patterns)

# rule WTS2_correlogram_τ50:
# 	output:
# 		# f"output/WTS2/{paramspace.wildcard_pattern}.png"
# 		expand('output/WTS2/ACFτ50/{params}_Acf.png', params = paramspace.wildcard_pattern),
# 		expand('output/WTS2/ACFτ50/{params}_pAcf.png', params = paramspace.wildcard_pattern)
# 	params:
# 		args1 = lambda w: w["SampleNumber"],
# 		args2 = lambda w: w["CellNumber"],
# 		args3 = lambda w: w["CellType"]
		
# 	benchmark:
# 		f'benchmarks/WTS2/ACFτ50/{paramspace.wildcard_pattern}.txt'
# 	conda:
# 		'envs/myenv_WTS2.yaml'
# 	resources:
# 		mem_gb=200
# 	log:
# 		f'logs/WTS2/ACFτ50/{paramspace.wildcard_pattern}.log'
# 	shell:
# 		'src/WTS2_correlogram.sh {params.args1} {params.args2} {params.args3} >& {log}'
# ###################################################

# # WTS2 correlogram τ500
# ###################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"

# # read sample_sheet
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')
# # remove sample
# SAMPLE_SHEET = SAMPLE_SHEET[~SAMPLE_SHEET['SampleNumber'].isin(['7','12','13','14'])]
# SAMPLE_SHEET  = SAMPLE_SHEET .reset_index(drop=True)


# # paramspace
# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
# 	input:
# 		expand('output/WTS2/ACFτ500/{params}_Acf.png', params = paramspace.instance_patterns),
# 		expand('output/WTS2/ACFτ500/{params}_pAcf.png', params = paramspace.instance_patterns)

# rule WTS2_correlogram_τ500:
# 	output:
# 		# f"output/WTS2/ACFτ500/{paramspace.wildcard_pattern}.png"
# 		expand('output/WTS2/ACFτ500/{params}_Acf.png', params = paramspace.wildcard_pattern),
# 		expand('output/WTS2/ACFτ500/{params}_pAcf.png', params = paramspace.wildcard_pattern)
# 	params:
# 		args1 = lambda w: w["SampleNumber"],
# 		args2 = lambda w: w["CellNumber"],
# 		args3 = lambda w: w["CellType"]
		
# 	benchmark:
# 		f'benchmarks/WTS2/ACFτ500/{paramspace.wildcard_pattern}.txt'
# 	conda:
# 		'envs/myenv_WTS2.yaml'
# 	resources:
# 		mem_gb=200
# 	log:
# 		f'logs/WTS2/ACFτ500/{paramspace.wildcard_pattern}.log'
# 	shell:
# 		'src/WTS2_correlogram.sh {params.args1} {params.args2} {params.args3} >& {log}'
# ###################################################

# # WTS2 correlogram stim τ50
# #####################################################################################
# import pandas as pd
# from snakemake.utils import Paramspace

# configfile: "config.yaml"

# # read sample_sheet
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')
# # remove sample
# SAMPLE_SHEET = SAMPLE_SHEET[~SAMPLE_SHEET['SampleNumber'].isin(['7','12','13','14'])]
# SAMPLE_SHEET  = SAMPLE_SHEET .reset_index(drop=True)

# # paramspace
# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
# 	input:
# 		expand('output/WTS2/ACFτ50_stim/ACF/{params}_BF.png', params = paramspace.instance_patterns),
# 		expand('output/WTS2/ACFτ50_stim/ACF/{params}_AF.png', params = paramspace.instance_patterns),
# 		expand('output/WTS2/ACFτ50_stim/pACF/{params}_BF.png', params = paramspace.instance_patterns),
# 		expand('output/WTS2/ACFτ50_stim/pACF/{params}_AF.png', params = paramspace.instance_patterns)

# rule WTS2_correlogram_τ50:
# 	output:
# 		expand('output/WTS2/ACFτ50_stim/ACF/{params}_BF.png', params = paramspace.wildcard_pattern),
# 		expand('output/WTS2/ACFτ50_stim/ACF/{params}_AF.png', params = paramspace.wildcard_pattern),
# 		expand('output/WTS2/ACFτ50_stim/pACF/{params}_BF.png', params = paramspace.wildcard_pattern),
# 		expand('output/WTS2/ACFτ50_stim/pACF/{params}_AF.png', params = paramspace.wildcard_pattern)
# 	params:
# 		args1 = lambda w: w["SampleNumber"],
# 		args2 = lambda w: w["CellNumber"],
# 		args3 = lambda w: w["CellType"]
		
# 	benchmark:
# 		f'benchmarks/WTS2/ACFτ50_stim/{paramspace.wildcard_pattern}.txt'
# 	conda:
# 		'envs/myenv_WTS2.yaml'
# 	resources:
# 		mem_gb=200
# 	log:
# 		f'logs/WTS2/ACFτ50_stim/{paramspace.wildcard_pattern}.log'
# 	shell:
# 		'src/WTS2_correlogram.sh {params.args1} {params.args2} {params.args3} >& {log}'
# #####################################################################################