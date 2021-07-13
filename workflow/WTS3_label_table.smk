# WTS3 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
# remove artifact
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')
# remove not annotation ASER, not periodic ASER →17sample
N_SAMPLES.remove('4')
# N_SAMPLES.remove('7')
N_SAMPLES.remove('9')
N_SAMPLES.remove('15')
N_SAMPLES.remove('16')
N_SAMPLES.remove('22')
N_SAMPLES.remove('26')
N_SAMPLES.remove('28')

# Distance Data
dist_data = ["SBD"]
# data time range
time_range =  ["all","stimAfter"]
# Clustering Evaluation Method
cls_eval = ["ARI"]
# cls_eval = ["purity", "ARI", "Fmeasure", "Entropy"]
# salt-stimulated associated cells
stim_label = ["label_acf", "label_cls"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/table/{label}/table_neuron.csv', 
            dist=dist_data,
            range=time_range,
            eval=cls_eval,
            label=stim_label
            ),
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/table/{label}/table_rmSensory.csv', 
            dist=dist_data,
            range=time_range,
            eval=cls_eval,
            label=stim_label
            ),
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/table/{label}/table_yshift.csv', 
            dist=dist_data,
            range=time_range,
            eval=cls_eval,
            label=stim_label
            )
rule label_table:
    # input:
    #     expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/SampleNumber_{N}/label_table.RData',N=N_SAMPLES),
    #     expand('output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_value.RData',N=N_SAMPLES)
    output:
        table_neuron = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/table/{label}/table_neuron.csv',
        table_rmSensory = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/table/{label}/table_rmSensory.csv',
        table_yshift = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/table/{label}/table_yshift.csv'
    params:
        label_dir = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/SampleNumber_',
        yshift_value_dir = 'output/WTS3/normalize_1/{range}/{dist}//SampleNumber_'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/{eval}/table/{label}/table.txt'
    conda:
        '../envs/myenv_WTS3_tsPlot.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/{eval}/table/{label}/table.log'
    shell:
        'src/WTS3_label_table.sh {params.label_dir} {params.yshift_value_dir} {wildcards.label} {output.table_neuron} {output.table_rmSensory} {output.table_yshift} >& {log}'
###################################################