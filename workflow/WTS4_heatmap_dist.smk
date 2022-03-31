# WTS4_heatmap_dist
##################################################
# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Dist_heatmap/SampleNumber_{N}.png',
            range=time_range,
            dist=dist_data,
            N=N_SAMPLES
            )
        
rule WTS4_heatmap_dist:
    input:
        distance = 'output/WTS4/normalize_1/{range}/{dist}/Distance/SampleNumber_{N}.RData'
    output:
        heatmap = 'output/WTS4/normalize_1/{range}/{dist}/Dist_heatmap/SampleNumber_{N}.png'
    params:
        merged_cls = 'output/WTS4/normalize_1/{range}/{dist}/MCMIHOOI/Merged_cls/k_Number_9.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Dist_heatmap/SampleNumber_{N}.txt'
    container:
        "docker://yamaken37/heatmap_dist:20220331"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Dist_heatmap/SampleNumber_{N}.log'
    shell:
        'src/WTS4_heatmap_dist.sh {input.distance} {params.merged_cls} {output.heatmap}>& {log}'