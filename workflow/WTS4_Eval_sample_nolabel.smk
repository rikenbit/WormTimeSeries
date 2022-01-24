# WTS4_Eval_sample_nolabel
###################################################

# No. of Clusters
# N_CLUSTERS = list(map(str, range(2, 21)))
N_CLUSTERS = ["3"]

# Distance Data
# dist_data = ["EUCL","SBD_abs"]
dist_data = ["SBD_abs"]

# data time range
time_range = ["stimAfter"]

# Evaluation Method
Evaluation_method = ["PseudoF","Connectivity","kNN"]

rule all:
    input: 
        expand('output/WTS4/normalize_1/{range}/{dist}/Eval_sample/{Eval}_nolabel/k_Number_{N_cls}.RData',
            range=time_range,
            dist=dist_data,
            N_cls=N_CLUSTERS,
            Eval=Evaluation_method
            )
        
rule WTS4_Eval_sample_nolabel:
    input:
        sample_cls =  'output/WTS4/normalize_1/{range}/{dist}/Cluster_sample/k_Number_{N_cls}/sample_cls.RData',
    output:
        eval_result = 'output/WTS4/normalize_1/{range}/{dist}/Eval_sample/{Eval}_nolabel/k_Number_{N_cls}.RData'
    params:
        input_sample  = 'output/WTS4/normalize_1/{range}/{dist}/Distance',
        input_path = 'data/normalize_1',
        stim_xlsx = 'data/stimulation/stimulation_timing.xlsx'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/Eval_sample/{Eval}_nolabel/k_Number_{N_cls}.txt'
    container:
        "docker://yamaken37/eval_sample:20220117"
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/Eval_sample/{Eval}_nolabel/k_Number_{N_cls}.log'
    shell:
        'src/WTS4_Eval_sample_nolabel.sh {params.input_sample} {input.sample_cls} {wildcards.range} {wildcards.N_cls} {wildcards.Eval} {params.input_path} {params.stim_xlsx} {output.eval_result} >& {log}'