# WTS4_MCMI_pairs
###################################################
# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["9"]

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["MCMIHOOI"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Cluster.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Neuron_type.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            ),
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Class.png',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Re_cls=ReClustering_method
            )
        
rule WTS4_MCMI_pairs:
    input:
        Mem_matrix = 'output/WTS4/normalize_1/{range}/{dist}/Membership/k_Number_{N_cls}.RData'
    output:
        Cluster = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Cluster.png',
        Neuron_type = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Neuron_type.png',
        Class = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/Class.png'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/MCMI_pairs.txt'
    container:
        "docker://docker_images"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/PairPlot/k_Number_{N_cls}/MCMI_pairs.log'
    shell:
        'src/WTS4_MCMI_pairs.sh {wildcards.Re_cls} {input.Mem_matrix} {output.Cluster} {output.Neuron_type} {output.Class}>& {log}'