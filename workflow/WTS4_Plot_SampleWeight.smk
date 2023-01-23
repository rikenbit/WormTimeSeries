# WTS4_Plot_SampleWeight
###################################################
# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# ReClustering method
ReClustering_method = ["MCMIHOOI"]
# ReClustering_method = ["CSPA","OINDSCAL","MCMIHOOI"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Plot_SampleWeight.svg',
            range=time_range,
            dist=dist_data,
            Re_cls=ReClustering_method
            )

rule WTS4_Plot_SampleWeight:
    output:
        'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Plot_SampleWeight.svg'
    params:
        weight_path  = 'output/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Merged_data',
        sample_path  = 'output/WTS4/normalize_1/{range}/{dist}/Distance'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Plot_SampleWeight.txt'
    container:
        # "docker://yamaken37/cluster_sample:20220119"
        "docker://yamaken37/ggplot_svg:20230118"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/{Re_cls}/Plot_SampleWeight.log'
    shell:
        'src/WTS4_Plot_SampleWeight.sh {params.weight_path} {params.sample_path} {output} >& {log}'