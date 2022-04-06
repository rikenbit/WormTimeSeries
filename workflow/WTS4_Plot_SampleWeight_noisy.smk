# WTS4_Plot_SampleWeight_noisy ラベルノイズ
###################################################
# NOISE SAMPLE PATERN
NOISE_TEST = ["n1_noise"]

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["MCMIHOOI"]

# NOISE TRY
N_TRY = list(map(str, range(1, 11)))
# N_TRY = ["1"]

# SAMPLE NUMBER
# N_SAMPLES = list(map(str, range(1, 29)))
# # remove artifact
# N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')
N_SAMPLES = ["1"]
rule all:
    input: 
        expand('output/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/{Re_cls}/Plot_SampleWeight.png',
            range=time_range,
            dist=dist_data,
            Re_cls=ReClustering_method,
            NOISE=NOISE_TEST,
            N=N_SAMPLES,
            TRY=N_TRY
            )
        
rule WTS4_Plot_SampleWeight_noisy:
    output:
        'output/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/{Re_cls}/Plot_SampleWeight.png'
    params:
        weight_path  = 'output/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/{Re_cls}/Merged_data',
        sample_path  = 'output/WTS4/{NOISE}_{TRY}/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/{Re_cls}/Plot_SampleWeight.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{NOISE}_{TRY}/NoiseSampleNumber_{N}/{range}/{dist}/{Re_cls}/Plot_SampleWeight.log'
    shell:
        'src/WTS4_Plot_SampleWeight.sh {params.weight_path} {params.sample_path} {output}>& {log}'

# # WTS4_Plot_SampleWeight_noisy
# ###################################################
# # NOISE SAMPLE PATERN
# # NOISE_TEST = ["normalize_1","n1_28sample","n1_24sample_add3","n1_24sample_add8","n1_24sample_add20","n1_24sample_add25"]
# # NOISE_TEST = ["n1_28sample"]
# # NOISE_TEST = ["normalize_1"]
# # NOISE_TEST = ["n1_28sample_trim20"]
# NOISE_TEST = ["n1_27sample_rm20"]

# # Distance Data
# # dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# # data time range
# time_range = ["stimAfter"]

# # ReClustering method
# ReClustering_method = ["MCMIHOOI"]

# rule all:
#     input: 
#         expand('output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Plot_SampleWeight.png',
#             range=time_range,
#             dist=dist_data,
#             Re_cls=ReClustering_method,
#             NOISE=NOISE_TEST
#             )
        
# rule WTS4_Plot_SampleWeight_noisy:
#     output:
#         'output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Plot_SampleWeight.png'
#     params:
#         weight_path  = 'output/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Merged_data',
#         sample_path  = 'output/WTS4/{NOISE}/{range}/{dist}/Distance'
#     benchmark:
#         'benchmarks/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Plot_SampleWeight.txt'
#     container:
#         "docker://yamaken37/cluster_sample:20220119"
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/{NOISE}/{range}/{dist}/{Re_cls}/Plot_SampleWeight.log'
#     shell:
#         'src/WTS4_Plot_SampleWeight.sh {params.weight_path} {params.sample_path} {output}>& {log}'