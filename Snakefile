import pandas as pd
from snakemake.utils import min_version
from snakemake.utils import Paramspace

# min_version("5.3.2")
configfile: "config.yaml"

N_SAMPLES = list(map(str, range(1, 16)))

# WTS1 sample_sheet
SAMPLE_SHEET = pd.read_csv(config['SAMPLE_SHEET'], dtype='string')
# test 5行
SAMPLE_SHEET = SAMPLE_SHEET[206:211]
SAMPLE_SHEET.columns = ['Samplenumber','Cellnumber','Celltype']

paramspace = Paramspace(SAMPLE_SHEET, filename_params=['Cellnumber', 'Celltype'], param_sep="")



rule all:
	input:
		expand('output/WTS1/{params}.png', params = paramspace.instance_patterns)

# WTS1 dataload
rule WTS1_load:
	input:
		expand('data/cleandata/{N}_ratio.csv', N=N_SAMPLES),
		'data/stimulation_timing.csv'
	output:
		expand('data/cleandata_mat/matrix_{N}.RData', N=N_SAMPLES),
		expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES),
		'data/WTS1_sample_sheet.csv'
	benchmark:
		'benchmarks/WTS1_load.txt'
	conda:
		'envs/myenv_WTS1.yaml'
	resources:
		mem_gb=200
	log:
		'logs/WTS1_load.log'
	shell:
		'src/WTS1_load.sh >& {log}'

# WTS1 1cell plot
rule WTS1_1cell:
	input:
		expand('data/cleandata_mat/matrix_{N}.RData', N=N_SAMPLES),
		expand('data/stimulation/stim_{N}.RData', N=N_SAMPLES)
	output:
		# 'output/WTS1/Samplenumber1/Cellnumber207_Celltypeg1P.png',
		# 'output/WTS1/Samplenumber1/Cellnumber208_Celltypeg2L.png',
		# 'output/WTS1/Samplenumber1/Cellnumber209_Celltypeg2R.png',
		# 'output/WTS1/Samplenumber2/Cellnumber1_CelltypeX103.png',
		# 'output/WTS1/Samplenumber2/Cellnumber2_CelltypeX109.png'
		# {input.all}
		# expand(f'output/WTS1/{paramspace.instance_patterns}.png', params = paramspace.instance_patterns)
		# expand(f'output/WTS1/{paramspace.wildcard_pattern}.png')
		# f'output/WTS1/{input.instance_patterns}.png'
		# expand(f'output/WTS1/{paramspace.wildcard_pattern}.png', Samplenumber = lambda w: w["Samplenumber"], Cellnumber = lambda w: w["Cellnumber"], Celltype = lambda w: w["Celltype"])
		# lambda w: 'output/WTS1/{w}.png'
		# f"output/WTS1/{paramspace.instance_patterns}.png"
		# f"output/WTS1/{paramspace.wildcard_pattern}.png"
		# f"expand('output/WTS1/{params}.png', params = paramspace.instance_patterns)"
		# expand('output/WTS1/{params}.png', params = paramspace.instance_patterns)
		f"output/WTS1/{paramspace.wildcard_pattern}.png"
	params:
		args1 = lambda w: w["Samplenumber"],
		args2 = lambda w: w["Cellnumber"],
		args3 = lambda w: w["Celltype"]
		
	benchmark:
		# expand('benchmarks/WTS1/{params}.txt',params = paramspace.wildcard_pattern)
		# 'benchmarks/WTS1/{wildcards.params}.txt'
		f'benchmarks/WTS1/{paramspace.wildcard_pattern}.txt'
	conda:
		'envs/myenv_WTS1.yaml'
	resources:
		mem_gb=200
	log:
		# expand('logs/WTS1/{params}.log',params = paramspace.wildcard_pattern)
		# 'logs/WTS1/{wildcards.params}.log'
		f'logs/WTS1/{paramspace.wildcard_pattern}.log'
	shell:
		'src/WTS1_1cell.sh {params.args1} {params.args2} {params.args3} >& {log}'