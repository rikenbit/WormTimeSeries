# WTS4 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# No. of Clusters とりあえずk=3？
N_CLUSTERS = list(map(str, range(3, 21)))

# Distance Data
dist_data = ["EUCL","SBD_abs"]
# data time range
time_range = ["all","stimAfter"]

rule all:
    input:
        expand('output/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{N_cls}_Clusters/Membership.RData', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range,
            N_cls=N_CLUSTERS
            )
        
rule Membership:
    input:
        dist_matrix = 'output/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{dist}.RData'
    output:
        Mem_matrix = 'output/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{N_cls}_Clusters/Membership.RData'
    benchmark:
        'benchmarks/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{N_cls}_Clusters/Membership.txt'
    conda:
        '../envs/myenv_WTS4_Membership.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS4/normalize_1/{range}/{dist}/SampleNumber_{N}/{N_cls}_Clusters/Membership.log'
    shell:
        'src/WTS4_Membership.sh {wildcards.N_cls} {input.dist_matrix} {output.Mem_matrix}>& {log}'