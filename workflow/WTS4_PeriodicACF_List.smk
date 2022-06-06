# WTS4_PeriodicACF_List
###################################################
# No. of observation
N_observation = ["5"]

# normalize pattern
normalize_pattern = ["normalize_1"]

rule all:
    input:
        expand('output/WTS4/{normalize_P}/PeriodicACF/periodic_acf_over{N_obs}.csv',
            N_obs=N_observation,
            normalize_P=normalize_pattern
            ),
        expand('output/WTS4/{normalize_P}/PeriodicACF/sample_over{N_obs}.csv',
            N_obs=N_observation,
            normalize_P=normalize_pattern
            )
        
rule WTS4_PeriodicACF_List:
    input:
        csv = 'output/WTS2/WTS2_PeriodicACF.csv'
    output:
        acf_list = 'output/WTS4/{normalize_P}/PeriodicACF/periodic_acf_over{N_obs}.csv',
        acf_sample = 'output/WTS4/{normalize_P}/PeriodicACF/sample_over{N_obs}.csv'
    benchmark:
        'benchmarks/WTS4/{normalize_P}/PeriodicACF/periodic_acf_over{N_obs}.txt'
    container:
        "docker://yamaken37/cluster_sample:20220119"
        # tidyverseが入っているコンテナ
    resources:
        mem_gb=200
    log:
        'logs/WTS4/{normalize_P}/PeriodicACF/periodic_acf_over{N_obs}.log'
    shell:
        'src/WTS4_PeriodicACF_List.sh {input.csv} {wildcards.N_obs} {output.acf_list} {output.acf_sample} >& {log}'