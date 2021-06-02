# WTS3 SBD 
###################################################
N_SAMPLES = list(map(str, range(1, 29)))
N_SAMPLES.remove('3')
N_SAMPLES.remove('8')
N_SAMPLES.remove('20')
N_SAMPLES.remove('25')

# Distance Data
dist_data = ["SBD"]

# Dimensionality Reduction Method
DR_Method = ["tsne","umap"]

#Clustering Evaluation Method
cls_eval = ["ARI", "Fmeasure"]

rule all:
    input:
        expand('output/WTS3/{dist}/normalize_1/all/{dr_method}/{eval}/plot/SampleNumber_{N}.png', 
            N=N_SAMPLES, 
            eval=cls_eval, 
            dr_method=DR_Method,
            dist=dist_data)

rule WTS3_plot:
    input:
        Neuron = 'output/data/normalize_1/ReadData_{N}.RData',
        stim = 'output/data/stimulation/stim_{N}.RData',
        mCherry = 'output/data/mCherry/mCherry_{N}.RData',
        Position = 'output/data/Position/Position_{N}.RData',
        tempdata = 'output/WTS3/{dist}/normalize_1/all/{dr_method}/{eval}/cls_tempdata/SampleNumber_{N}.RData'
        
    output:
        png = 'output/WTS3/{dist}/normalize_1/all/{dr_method}/{eval}/plot/SampleNumber_{N}.png'
    benchmark:
        'benchmarks/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.txt'
    conda:
        '../envs/myenv_WTS3.yaml'
    resources:
        mem_gb=200
    log:
        'logs/WTS3/SBD/normalize_1/all/SampleNumber_{N}/SBD.log'
    shell:
        'src/WTS3_SBD.sh {wildcards.N} {input.Neuron} {input.stim} {input.mCherry} {input.Position} {input.tempdata} {output.png} {wildcards.eval} {wildcards.dr_method} >& {log}'
###################################################