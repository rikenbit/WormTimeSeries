# WTS1 dataload
# ###################################################
# N_SAMPLES = list(map(str, range(1, 16)))
# 
# rule all:
# 	input:
# 		expand('data/cleandata_mat/matrix_{N}.RData', N=N_SAMPLES),
#         expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES),
#         'data/WTS1_sample_sheet.csv'
# 
# rule WTS1_load:
#     input:
#         expand('data/cleandata/{N}_ratio.csv', N=N_SAMPLES),
#         'data/stimulation_timing.csv'
#     output:
#         expand('data/cleandata_mat/matrix_{N}.RData', N=N_SAMPLES),
#         expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES),
#         'data/WTS1_sample_sheet.csv'
#     benchmark:
#         'benchmarks/WTS1/WTS1_load.txt'
#     conda:
#         'envs/myenv_WTS1.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS1/WTS1_load.log'
#     shell:
#         'src/WTS1_load.sh >& {log}'
# ###################################################

# # WTS1 1cell plot
# ###################################################
# import pandas as pd
# # from snakemake.utils import min_version
# from snakemake.utils import Paramspace


# # min_version("5.3.2")
# configfile: "config.yaml"

# N_SAMPLES = list(map(str, range(1, 16)))

# # read sample_sheet
# SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')

# # paramspace
# paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

# rule all:
# 	input:
# 		expand('output/WTS1/{params}.png', params = paramspace.instance_patterns)

# rule WTS1_1cell:
# 	# input:
# 	# 	expand('data/cleandata_mat/matrix_{N}.RData', N=N_SAMPLES),
#  #        expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES),
#  #        'data/WTS1_sample_sheet.csv'
# 	output:
# 		f"output/WTS1/{paramspace.wildcard_pattern}.png"
# 	params:
# 		args1 = lambda w: w["SampleNumber"],
# 		args2 = lambda w: w["CellNumber"],
# 		args3 = lambda w: w["CellType"]
		
# 	benchmark:
# 		f'benchmarks/WTS1/{paramspace.wildcard_pattern}.txt'
# 	conda:
# 		'envs/myenv_WTS1.yaml'
# 	resources:
# 		mem_gb=200
# 	log:
# 		f'logs/WTS1/{paramspace.wildcard_pattern}.log'
# 	shell:
# 		'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} >& {log}'
# ###################################################

# WTS2 correlogram
###################################################
import pandas as pd
from snakemake.utils import Paramspace

configfile: "config.yaml"

# read sample_sheet
SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')
# remove sample
SAMPLE_SHEET = SAMPLE_SHEET[~SAMPLE_SHEET['SampleNumber'].isin(['7','12','13','14'])]
SAMPLE_SHEET  = SAMPLE_SHEET .reset_index(drop=True)

# test
###################################################
SAMPLE_SHEET = SAMPLE_SHEET[206:211]
###################################################

# paramspace
paramspace = Paramspace(SAMPLE_SHEET, filename_params=['CellNumber', 'CellType'], param_sep="_")

rule all:
	input:
		expand('output/WTS2/{params}_Acf.png', params = paramspace.instance_patterns),
		expand('output/WTS2/{params}_pAcf.png', params = paramspace.instance_patterns)

rule WTS2_correlogram:
	output:
		# f"output/WTS2/{paramspace.wildcard_pattern}.png"
		expand('output/WTS2/{params}_Acf.png', params = paramspace.wildcard_pattern),
		expand('output/WTS2/{params}_pAcf.png', params = paramspace.wildcard_pattern)
	params:
		args1 = lambda w: w["SampleNumber"],
		args2 = lambda w: w["CellNumber"],
		args3 = lambda w: w["CellType"]
		
	benchmark:
		f'benchmarks/WTS2/{paramspace.wildcard_pattern}.txt'
	conda:
		'envs/myenv_WTS2.yaml'
	resources:
		mem_gb=200
	log:
		f'logs/WTS2/{paramspace.wildcard_pattern}.log'
	shell:
		'src/WTS2_correlogram.sh {params.args1} {params.args2} {params.args3} >& {log}'
###################################################