# # WTS4 SBD 
# ###################################################
# N_SAMPLES = list(map(str, range(1, 29)))

# # remove artifact
# N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')


# # Distance Data
# dist_data = ["EUCL","SBD_abs"]
# # data time range
# time_range = ["all","stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData', 
#             N=N_SAMPLES,
#             dist=dist_data,
#             range=time_range
#             )
        
# rule DistMatrix:
#     input:
#         RData = 'data/normalize_1/ReadData_{N}.RData'
#     output:
#         dist_matrix = 'output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData'
#     params:
#         stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
#     benchmark:
#         'benchmarks/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.txt'
#     conda:
#         '../envs/myenv_WTS4_{dist}.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.log'
#     shell:
#         'src/WTS4_{wildcards.dist}.sh {wildcards.N} {input.RData} {wildcards.range} {params.stim_xlsx} {output.dist_matrix}>& {log}'

# # WTS4 SBD 
# ###################################################
# N_SAMPLES = ["3","8","20","25"]

# # NOISE_TEST = ["n1_28sample","normalize_1"]
# NOISE_TEST = ["n1_28sample"]

# # Distance Data
# dist_data = ["EUCL","SBD_abs"]
# # data time range
# time_range = ["all","stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/{NOISE}/{range}/{dist}/Distance/SampleNumber_{N}.RData', 
#             N=N_SAMPLES,
#             dist=dist_data,
#             range=time_range,
#             NOISE=NOISE_TEST
#             )
        
# rule DistMatrix:
#     input:
#         RData = 'data/{NOISE}/ReadData_{N}.RData'
#     output:
#         dist_matrix = 'output/WTS4/{NOISE}/{range}/{dist}/Distance/SampleNumber_{N}.RData'
#     params:
#         stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
#     benchmark:
#         'benchmarks/WTS4/{NOISE}/{range}/{dist}/Distance/SampleNumber_{N}.txt'
#     conda:
#         '../envs/myenv_WTS4_{dist}.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/{NOISE}/{range}/{dist}/Distance/SampleNumber_{N}.log'
#     shell:
#         'src/WTS4_{wildcards.dist}.sh {wildcards.N} {input.RData} {wildcards.range} {params.stim_xlsx} {output.dist_matrix}>& {log}'

# # WTS4 SBD n1_28sample_trim20
# ###################################################

# N_SAMPLES = ["20"]

# # Distance Data
# dist_data = ["SBD_abs"]
# # data time range
# time_range = ["stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/n1_28sample_trim20/{range}/{dist}/Distance/SampleNumber_{N}.RData', 
#             N=N_SAMPLES,
#             dist=dist_data,
#             range=time_range
#             )
        
# rule WTS4_DistMatrix:
#     input:
#         RData = 'data/normalize_1/ReadData_{N}.RData'
#     output:
#         dist_matrix = 'output/WTS4/n1_28sample_trim20/{range}/{dist}/Distance/SampleNumber_{N}.RData'
#     params:
#         stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
#     benchmark:
#         'benchmarks/WTS4/n1_28sample_trim20/{range}/{dist}/Distance/SampleNumber_{N}.txt'
#     conda:
#         '../envs/myenv_WTS4_{dist}.yaml'
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/n1_28sample_trim20/{range}/{dist}/Distance/SampleNumber_{N}.log'
#     shell:
#         'src/WTS4_{wildcards.dist}.sh {wildcards.N} {input.RData} {wildcards.range} {params.stim_xlsx} {output.dist_matrix}>& {log}'

# WTS4 SBD n1_noise
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
# N_SAMPLES = ["2"]

# data group
NOISE_TEST = ["n1_noise", "n1_noise_TF"]

N_TRY = list(map(str, range(1, 11)))
# N_TRY = ["3"]

# Distance Data
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

rule all:
    input:
        expand('output/WTS4/{NOISE}_{TRY}/{range}/{dist}/Distance/SampleNumber_{N}.RData', 
            N=N_SAMPLES,
            NOISE=NOISE_TEST,
            TRY=N_TRY,
            dist=dist_data,
            range=time_range
            )
        
rule WTS4_DistMatrix:
    input:
        RData = 'data/{NOISE}_{TRY}/ReadData_{N}.RData'
    output:
        dist_matrix = 'output/WTS4/{NOISE}_{TRY}/{range}/{dist}/Distance/SampleNumber_{N}.RData'
    params:
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS4/{NOISE}_{TRY}/{range}/{dist}/Distance/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS4_{dist}.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{NOISE}_{TRY}/{range}/{dist}/Distance/SampleNumber_{N}.log'
    shell:
        'src/WTS4_{wildcards.dist}.sh {wildcards.N} {input.RData} {wildcards.range} {params.stim_xlsx} {output.dist_matrix}>& {log}'