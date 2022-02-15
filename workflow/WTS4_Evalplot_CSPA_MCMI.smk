# WTS4_Evalplot_CSPA_MCMI
###################################################
# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# Evaluation Method
Evaluation_method = ["PseudoF", "Connectivity", "kNN", "ARI_behavior", "purity_behavior", "Fmeasure_behavior", "Entropy_behavior", "NMI_behavior", "AMI_behavior","Silhouette"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/Evalplot_CSPA_MCMI/{Eval}.png',
            range=time_range,
            dist=dist_data,
            Eval=Evaluation_method
            )
        
rule WTS4_Evalplot_CSPA_MCMI:
    output:
        'output/WTS4/normalize_1/{range}/{dist}/Evalplot_CSPA_MCMI/{Eval}.png'
    params:
        input_value  = 'output/WTS4/normalize_1/{range}/{dist}'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Evalplot_CSPA_MCMI/{Eval}.txt'
    container:
        "docker://yamaken37/dimreduc_mcmi:20211224"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Evalplot_CSPA_MCMI/{Eval}.log'
    shell:
        'src/WTS4_Evalplot_CSPA_MCMI.sh {params.input_value} {wildcards.Eval} {output}>& {log}'