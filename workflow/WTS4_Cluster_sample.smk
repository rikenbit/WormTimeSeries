# WTS4_Cluster_sample
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# No. of Clusters
N_CLUSTERS = list(map(str, range(2, 21)))
# N_CLUSTERS = ["3"]

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS
            )
        
rule WTS4_Cluster_sample:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData',
            dist=dist_data,
            range=time_range,
            N=N_SAMPLES
            )
    params:
        dist_matrix_path = 'output/WTS4/normalize_1/{range}/{dist}/Distance'
    output:
        'output/WTS4/normalize_1/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.txt'
    container:
        "docker:docker_images"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.log'
    shell:
        'src/WTS4_Cluster_sample.sh {params.dist_matrix_path} {output} {wildcards.N_cls} >& {log}'