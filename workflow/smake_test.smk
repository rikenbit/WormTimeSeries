# smake_test
###################################################
# No. of Clusters
# test_number = list(map(str, range(2, 21)))
test_number = ["3"]


rule all:
    input:
        expand('data/smake_test/{TEST_N}_df.csv', TEST_N=test_number)

rule smake_test:
    input:
        'data/smake_test/input_number.csv'
    output:
        'data/smake_test/{TEST_N}_df.csv'
    benchmark:
        'benchmarks/smake_test/{TEST_N}_df.txt'
    container:
        "docker://yamaken37/ggplot_svg:20230118"
    resources:
        mem_gb=200
    log:
        'logs/smake_test/{TEST_N}_df.log'
    shell:
        'src/smake_test.sh {wildcards.TEST_N} {input} {output} >& {log}'