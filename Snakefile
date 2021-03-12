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

# # WTS1 plot
# ###################################################
# import pandas as pd
# # from snakemake.utils import min_version
# from snakemake.utils import Paramspace

# # min_version("5.3.2")
# configfile: "config.yaml"

# # read sample_sheet
# DATA_DIR_CFP = config['DATA_DIR_CFP']
# SAMPLE_SHEET_CFP = pd.read_csv(config['SAMPLE_SHEET_CFP'], dtype='string')
# paramspace_CFP = Paramspace(SAMPLE_SHEET_CFP, filename_params=['CellNumber', 'CellType'], param_sep="_")
# DATA_DIR_YFP = config['DATA_DIR_YFP']
# SAMPLE_SHEET_YFP = pd.read_csv(config['SAMPLE_SHEET_YFP'], dtype='string')
# paramspace_YFP = Paramspace(SAMPLE_SHEET_YFP, filename_params=['CellNumber', 'CellType'], param_sep="_")
# DATA_DIR_normalize_1 = config['DATA_DIR_normalize_1']
# SAMPLE_SHEET_normalize_1 = pd.read_csv(config['SAMPLE_SHEET_normalize_1'], dtype='string')
# paramspace_normalize_1 = Paramspace(SAMPLE_SHEET_normalize_1, filename_params=['CellNumber', 'CellType'], param_sep="_")
# DATA_DIR_normalize_2 = config['DATA_DIR_normalize_2']
# SAMPLE_SHEET_normalize_2 = pd.read_csv(config['SAMPLE_SHEET_normalize_2'], dtype='string')
# paramspace_normalize_2 = Paramspace(SAMPLE_SHEET_normalize_2, filename_params=['CellNumber', 'CellType'], param_sep="_")
# DATA_DIR_normalize_3 = config['DATA_DIR_normalize_3']
# SAMPLE_SHEET_normalize_3 = pd.read_csv(config['SAMPLE_SHEET_normalize_3'], dtype='string')
# paramspace_normalize_3 = Paramspace(SAMPLE_SHEET_normalize_3, filename_params=['CellNumber', 'CellType'], param_sep="_")
# DATA_DIR_normalize_4 = config['DATA_DIR_normalize_4']
# SAMPLE_SHEET_normalize_4 = pd.read_csv(config['SAMPLE_SHEET_normalize_4'], dtype='string')
# paramspace_normalize_4 = Paramspace(SAMPLE_SHEET_normalize_4, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
#     input:
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace_CFP.instance_patterns, DATA_DIR = DATA_DIR_CFP),
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace_YFP.instance_patterns, DATA_DIR = DATA_DIR_YFP),
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace_normalize_1.instance_patterns, DATA_DIR = DATA_DIR_normalize_1),
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace_normalize_2.instance_patterns, DATA_DIR = DATA_DIR_normalize_2),
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace_normalize_3.instance_patterns, DATA_DIR = DATA_DIR_normalize_3),
#         expand('output/WTS1/plot/{DATA_DIR}/{params}.png', params = paramspace_normalize_4.instance_patterns, DATA_DIR = DATA_DIR_normalize_4)


# rule WTS1_plot_CFP:
#     output:
#         f"output/WTS1/plot/{DATA_DIR_CFP}/{paramspace_CFP.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = DATA_DIR_CFP
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR_CFP}/{paramspace_CFP.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR_CFP}/{paramspace_CFP.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}'
# rule WTS1_plot_YFP:
#     output:
#         f"output/WTS1/plot/{DATA_DIR_YFP}/{paramspace_YFP.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = DATA_DIR_YFP
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR_YFP}/{paramspace_YFP.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR_YFP}/{paramspace_YFP.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}'
# rule WTS1_plot_normalize_1:
#     output:
#         f"output/WTS1/plot/{DATA_DIR_normalize_1}/{paramspace_normalize_1.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = DATA_DIR_normalize_1
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR_normalize_1}/{paramspace_normalize_1.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR_normalize_1}/{paramspace_normalize_1.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}'
# rule WTS1_plot_normalize_2:
#     output:
#         f"output/WTS1/plot/{DATA_DIR_normalize_2}/{paramspace_normalize_2.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = DATA_DIR_normalize_2
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR_normalize_2}/{paramspace_normalize_2.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR_normalize_2}/{paramspace_normalize_2.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}'
# rule WTS1_plot_normalize_3:
#     output:
#         f"output/WTS1/plot/{DATA_DIR_normalize_3}/{paramspace_normalize_3.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = DATA_DIR_normalize_3
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR_normalize_3}/{paramspace_normalize_3.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR_normalize_3}/{paramspace_normalize_3.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}'
# rule WTS1_plot_normalize_4:
#     output:
#         f"output/WTS1/plot/{DATA_DIR_normalize_4}/{paramspace_normalize_4.wildcard_pattern}.png"
#     params:
#         args1 = lambda w: w["SampleNumber"],
#         args2 = lambda w: w["CellNumber"],
#         args3 = lambda w: w["CellType"],
#         args4 = DATA_DIR_normalize_4
#     benchmark:
#         f'benchmarks/WTS1/plot/{DATA_DIR_normalize_4}/{paramspace_normalize_4.wildcard_pattern}.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         f'logs/WTS1/plot/{DATA_DIR_normalize_4}/{paramspace_normalize_4.wildcard_pattern}.log'
#     shell:
#         'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}'
# ###################################################
# test
###################################################
rule all:
    input:
        "output/WTS1/plot/normalize_4/SampleNumber_25/CellNumber_65_CellType_ASKR.png"

rule WTS1_plot_test:
    output:
        "output/WTS1/plot/normalize_4/SampleNumber_25/CellNumber_65_CellType_ASKR.png"
    params:
        args1 = "25",
        args2 = "65",
        args3 = "ASKR",
        args4 = "normalize_4"
    benchmark:
        f'benchmarks/WTS1/plot/normalize_4/SampleNumber_25/CellNumber_65_CellType_ASKR.txt'
    conda:
        'envs/myenv_WTS1.yaml'
    resources:
        mem_gb=200
    log:
        f'logs/WTS1/plot/normalize_4/SampleNumber_25/CellNumber_65_CellType_ASKR.log'
    shell:
        'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} {params.args4} >& {log}' 
###################################################

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