import pandas as pd
from snakemake.utils import min_version
# from snakemake.utils import Paramspace

min_version("5.3.2")
configfile: "config.yaml"

N_SAMPLES = list(map(str, range(1, 16)))

# WTS1 sample_sheet
SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')
# test 5行
SAMPLE_SHEET = SAMPLE_SHEET[206:211]
print(SAMPLE_SHEET)
# print(args1)
# print(args2)
# print(args3)
SAMPLE_NAME = Paramspace(SAMPLE_SHEET, filename_params=['Sample.number', 'Cell.number', 'Cell.type']) 
print(SAMPLE_NAME)


# rule all:
# 	input:
# 		expand('output/WTS1/celegans{args1}/{args1}_{args2}_{args3}.png',
# 				args1 = SAMPLE_SHEET['Sample.number'],
# 				args2 = SAMPLE_SHEET['Cell.number'],
# 				args3 = SAMPLE_SHEET['Cell.type']
# 				)

# # WTS1 dataload
# rule WTS1_load:
# 	input:
# 		expand('data/cleandata/{N}_ratio.csv', N=N_SAMPLES),
# 		'data/stimulation_timing.csv'
# 	output:
# 		expand('data/cleandata_mat/matrix_{N}.RData', N=N_SAMPLES),
# 		expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES),
# 		'data/WTS1_sample_sheet.csv'
# 	benchmark:
# 		'benchmarks/WTS1_load.txt'
# 	conda:
# 		'envs/myenv_WTS1.yaml'
# 	resources:
# 		mem_gb=200
# 	log:
# 		'logs/WTS1_load.log'
# 	shell:
# 		'src/WTS1_load.sh >& {log}'

# # WTS1 1cell plot
# rule WTS1_1cell:
# 	input:
# 		'data/cleandata_mat/matrix_{args1}.RData',
# 		'data/stimulation/stim_{args1}.RData'
# 	output:
# 		touch('output/WTS1/celegans{args1}/{args1}_{args2}_{args3}.png')
# 	benchmark:
# 		'benchmarks/WTS1_{args1}_{args2}_{args3}.txt'
# 	conda:
# 		'envs/myenv_WTS1.yaml'
# 	resources:
# 		mem_gb=200
# 	log:
# 		'logs/WTS1_{args1}_{args2}_{args3}.log'
# 	shell:
# 		'src/WTS1_1cell.sh {wildcards.args1} {wildcards.args2} {wildcards.args3} >& {log}'
