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
time_range = ["all"]

# Clustering Evaluation Method
cls_eval = ["ARI"]
# cls_eval = ["purity", "ARI", "Fmeasure", "Entropy"]

# Dimensionality Reduction Method
DimReduc = ["tsne"]
# DimReduc = ["tsne","umap"]

# salt-stimulated associated cells
stim_label = ["label_acf", "label_cls"]

rule all:
    input:
        expand('output/WTS3/normalize_1/{range}/{dist}/{eval}/tsPlot/{label}/SampleNumber_{N}.png', 
            N=N_SAMPLES,
            dist=dist_data,
            range=time_range,
            eval=cls_eval,
            label=stim_label
            )
        
rule tsPlot:
    input:
        input_n = 'data/normalize_1/ReadData_{N}.RData',
        input_stim = 'data/stimulation/stim_{N}.RData',
        input_mCherry = 'data/mCherry/mCherry_{N}.RData',
        input_Position = 'data/Position/Position_{N}.RData',
        yshift = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift.RData',
        yshift_value = 'output/WTS3/normalize_1/{range}/{dist}/SampleNumber_{N}/yshift_value.RData',
        label_table = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/SampleNumber_{N}/label_table.RData'
    output:
        tsPlot = 'output/WTS3/normalize_1/{range}/{dist}/{eval}/tsPlot/{label}/SampleNumber_{N}.png'
    params:
        args_shift = 'ASER'
    benchmark:
        'benchmarks/WTS3/normalize_1/{range}/{dist}/{eval}/tsPlot/{label}/SampleNumber_{N}.txt'
    conda:
        '../envs/myenv_WTS3_tsPlot.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/normalize_1/{range}/{dist}/{eval}/tsPlot/{label}/SampleNumber_{N}.log'
    shell:
        'src/WTS3_tsPlot.sh {wildcards.N} {input.input_n} {input.input_stim} {input.input_mCherry} {input.input_Position} {input.yshift} {input.yshift_value} {input.label_table} {wildcards.label} {params.args_shift} {output.tsPlot} >& {log}'
###################################################