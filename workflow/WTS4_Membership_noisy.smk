# WTS4_Membership n1_27sample_rm20
###################################################
# 28SAMPLES
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('20')
# N_SAMPLES = ["3"]

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# Distance Data
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

rule all:
    input:
        expand('output/WTS4/n1_27sample_rm20/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
            dist=dist_data,
            range=time_range,
            N_cls=N_CLUSTERS
            )
        
rule WTS4_n1_27sample_rm20:
    input:
        expand('output/WTS4/n1_27sample_rm20/{range}/{dist}/Distance/SampleNumber_{N}.RData',
            dist=dist_data,
            range=time_range,
            N=N_SAMPLES
            )
    output:
        Mem_matrix = 'output/WTS4/n1_27sample_rm20/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
    params:
        dist_matrix_path = 'output/WTS4/n1_27sample_rm20/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/n1_27sample_rm20/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
    # conda:
    #     '../envs/myenv_WTS4_Membership.yaml'
    container:
        "docker://yamaken37/cluster_sample:20220119"
        # tidyverseが入っているコンテナ
    resources:
        mem_gb=200
    log:
        'logs/WTS4/n1_27sample_rm20/{range}/{dist}/Membership/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Membership.sh {wildcards.N_cls} {params.dist_matrix_path} {output.Mem_matrix}>& {log}'
# # WTS4_Membership n1_28sample_trim20
# ###################################################
# # 28SAMPLES
# N_SAMPLES = list(map(str, range(1, 29)))
# # N_SAMPLES = ["3"]

# # No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
# # N_CLUSTERS = ["3"]

# # Distance Data
# dist_data = ["SBD_abs"]

# # data time range
# time_range = ["stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/n1_28sample_trim20/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
#             dist=dist_data,
#             range=time_range,
#             N_cls=N_CLUSTERS
#             )
        
# rule WTS4_n1_28sample_trim20:
#     input:
#         expand('output/WTS4/n1_28sample_trim20/{range}/{dist}/Distance/SampleNumber_{N}.RData',
#             dist=dist_data,
#             range=time_range,
#             N=N_SAMPLES
#             )
#     output:
#         Mem_matrix = 'output/WTS4/n1_28sample_trim20/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
#     params:
#         dist_matrix_path = 'output/WTS4/n1_28sample_trim20/{range}/{dist}/Distance'
#     benchmark:
#         'benchmarks/WTS4/n1_28sample_trim20/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
#     # conda:
#     #     '../envs/myenv_WTS4_Membership.yaml'
#     container:
#         "docker://yamaken37/cluster_sample:20220119"
#         # tidyverseが入っているコンテナ
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/n1_28sample_trim20/{range}/{dist}/Membership/k_Number_{N_cls}.log'
#     shell:
#         'src/WTS4_Membership.sh {wildcards.N_cls} {params.dist_matrix_path} {output.Mem_matrix}>& {log}'

# # WTS4_Membership n1_28sample
# ###################################################
# # NOISE SAMPLE PATERN
# # NOISE_TEST = ["n1_28sample","n1_24sample_add3","n1_24sample_add8","n1_24sample_add20","n1_24sample_add25"]

# # SAMPLES
# # 28SAMPLES
# N_SAMPLES = list(map(str, range(1, 29)))
# # # remove artifact
# # N_SAMPLES.remove('3')
# # N_SAMPLES.remove('8')
# # N_SAMPLES.remove('20')
# # N_SAMPLES.remove('25')

# # No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))

# # Distance Data
# dist_data = ["EUCL","SBD_abs"]

# # data time range
# time_range = ["all","stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/n1_28sample/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
#             dist=dist_data,
#             range=time_range,
#             N_cls=N_CLUSTERS
#             )
        
# rule WTS4_n1_28sample:
#     input:
#         expand('output/WTS4/n1_28sample/{range}/{dist}/Distance/SampleNumber_{N}.RData',
#             dist=dist_data,
#             range=time_range,
#             N=N_SAMPLES
#             )
#     output:
#         Mem_matrix = 'output/WTS4/n1_28sample/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
#     params:
#         dist_matrix_path = 'output/WTS4/n1_28sample/{range}/{dist}/Distance'
#     benchmark:
#         'benchmarks/WTS4/n1_28sample/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
#     # conda:
#     #     '../envs/myenv_WTS4_Membership.yaml'
#     container:
#         "docker://yamaken37/cluster_sample:20220119"
#         # tidyverseが入っているコンテナ
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/n1_28sample/{range}/{dist}/Membership/k_Number_{N_cls}.log'
#     shell:
#         'src/WTS4_Membership.sh {wildcards.N_cls} {params.dist_matrix_path} {output.Mem_matrix}>& {log}'

# # WTS4_Membership n1_24sample_add3
# ###################################################
# # SAMPLES
# # 28SAMPLES
# N_SAMPLES = list(map(str, range(1, 29)))
# # # remove artifact
# # N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')

# # No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))

# # Distance Data
# dist_data = ["EUCL","SBD_abs"]

# # data time range
# time_range = ["all","stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/n1_24sample_add3/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
#             dist=dist_data,
#             range=time_range,
#             N_cls=N_CLUSTERS
#             )
        
# rule WTS4_n1_24sample_add3:
#     input:
#         expand('output/WTS4/n1_24sample_add3/{range}/{dist}/Distance/SampleNumber_{N}.RData',
#             dist=dist_data,
#             range=time_range,
#             N=N_SAMPLES
#             )
#     output:
#         Mem_matrix = 'output/WTS4/n1_24sample_add3/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
#     params:
#         dist_matrix_path = 'output/WTS4/n1_24sample_add3/{range}/{dist}/Distance'
#     benchmark:
#         'benchmarks/WTS4/n1_24sample_add3/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
#     # conda:
#     #     '../envs/myenv_WTS4_Membership.yaml'
#     container:
#         "docker://yamaken37/cluster_sample:20220119"
#         # tidyverseが入っているコンテナ
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/n1_24sample_add3/{range}/{dist}/Membership/k_Number_{N_cls}.log'
#     shell:
#         'src/WTS4_Membership.sh {wildcards.N_cls} {params.dist_matrix_path} {output.Mem_matrix}>& {log}'

# # WTS4_Membership n1_24sample_add8
# ###################################################
# # SAMPLES
# # 28SAMPLES
# N_SAMPLES = list(map(str, range(1, 29)))
# # # remove artifact
# N_SAMPLES.remove('3')
# # N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')

# # No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))

# # Distance Data
# dist_data = ["EUCL","SBD_abs"]

# # data time range
# time_range = ["all","stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/n1_24sample_add8/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
#             dist=dist_data,
#             range=time_range,
#             N_cls=N_CLUSTERS
#             )
        
# rule WTS4_n1_24sample_add8:
#     input:
#         expand('output/WTS4/n1_24sample_add8/{range}/{dist}/Distance/SampleNumber_{N}.RData',
#             dist=dist_data,
#             range=time_range,
#             N=N_SAMPLES
#             )
#     output:
#         Mem_matrix = 'output/WTS4/n1_24sample_add8/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
#     params:
#         dist_matrix_path = 'output/WTS4/n1_24sample_add8/{range}/{dist}/Distance'
#     benchmark:
#         'benchmarks/WTS4/n1_24sample_add8/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
#     # conda:
#     #     '../envs/myenv_WTS4_Membership.yaml'
#     container:
#         "docker://yamaken37/cluster_sample:20220119"
#         # tidyverseが入っているコンテナ
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/n1_24sample_add8/{range}/{dist}/Membership/k_Number_{N_cls}.log'
#     shell:
#         'src/WTS4_Membership.sh {wildcards.N_cls} {params.dist_matrix_path} {output.Mem_matrix}>& {log}'

# # WTS4_Membership n1_24sample_add20
# ###################################################
# # SAMPLES
# # 28SAMPLES
# N_SAMPLES = list(map(str, range(1, 29)))
# # # remove artifact
# N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# # N_SAMPLES.remove('20')
# N_SAMPLES.remove('25')

# # No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))

# # Distance Data
# dist_data = ["EUCL","SBD_abs"]

# # data time range
# time_range = ["all","stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/n1_24sample_add20/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
#             dist=dist_data,
#             range=time_range,
#             N_cls=N_CLUSTERS
#             )
        
# rule WTS4_n1_24sample_add20:
#     input:
#         expand('output/WTS4/n1_24sample_add20/{range}/{dist}/Distance/SampleNumber_{N}.RData',
#             dist=dist_data,
#             range=time_range,
#             N=N_SAMPLES
#             )
#     output:
#         Mem_matrix = 'output/WTS4/n1_24sample_add20/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
#     params:
#         dist_matrix_path = 'output/WTS4/n1_24sample_add20/{range}/{dist}/Distance'
#     benchmark:
#         'benchmarks/WTS4/n1_24sample_add20/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
#     # conda:
#     #     '../envs/myenv_WTS4_Membership.yaml'
#     container:
#         "docker://yamaken37/cluster_sample:20220119"
#         # tidyverseが入っているコンテナ
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/n1_24sample_add20/{range}/{dist}/Membership/k_Number_{N_cls}.log'
#     shell:
#         'src/WTS4_Membership.sh {wildcards.N_cls} {params.dist_matrix_path} {output.Mem_matrix}>& {log}'

# # WTS4_Membership n1_24sample_add25
# ###################################################
# # SAMPLES
# # 28SAMPLES
# N_SAMPLES = list(map(str, range(1, 29)))
# # # remove artifact
# N_SAMPLES.remove('3')
# N_SAMPLES.remove('8')
# N_SAMPLES.remove('20')
# # N_SAMPLES.remove('25')

# # No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))

# # Distance Data
# dist_data = ["EUCL","SBD_abs"]

# # data time range
# time_range = ["all","stimAfter"]

# rule all:
#     input:
#         expand('output/WTS4/n1_24sample_add25/{range}/{dist}/Membership/k_Number_{N_cls}.RData', 
#             dist=dist_data,
#             range=time_range,
#             N_cls=N_CLUSTERS
#             )
        
# rule WTS4_n1_24sample_add25:
#     input:
#         expand('output/WTS4/n1_24sample_add25/{range}/{dist}/Distance/SampleNumber_{N}.RData',
#             dist=dist_data,
#             range=time_range,
#             N=N_SAMPLES
#             )
#     output:
#         Mem_matrix = 'output/WTS4/n1_24sample_add25/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
#     params:
#         dist_matrix_path = 'output/WTS4/n1_24sample_add25/{range}/{dist}/Distance'
#     benchmark:
#         'benchmarks/WTS4/n1_24sample_add25/{range}/{dist}/Membership/k_Number_{N_cls}.txt'
#     # conda:
#     #     '../envs/myenv_WTS4_Membership.yaml'
#     container:
#         "docker://yamaken37/cluster_sample:20220119"
#         # tidyverseが入っているコンテナ
#     resources:
#         mem_gb=200
#     log:
#         'logs/WTS4/n1_24sample_add25/{range}/{dist}/Membership/k_Number_{N_cls}.log'
#     shell:
#         'src/WTS4_Membership.sh {wildcards.N_cls} {params.dist_matrix_path} {output.Mem_matrix}>& {log}'