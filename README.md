# WormTimeSeries
- Since the neural activity patterns showed little correlation when viewed at the same time, we investigated whether there was any correlation with neural activity shifted by time τ.
- There are 7 analysis.
1. Time series neural activity plot (each cell)
2. Correlogram (each cell, each c.elegans)
3. Time series dimensionality reduction and clustering (each c.elegans)
4. Spectrogram (each cell)
5. Time-delay embedding (each cell)
6. Eigenvalue distribution of the history matrix of all cells (each c.elegans)
7. Recurrence plot (each cell)


Installation
======
```bash
git clone https://github.com/rikenbit/WormTimeSeries/
cd WormTimeSeries
snakemake -j 1
```
# Prerequisites
- Snakemake (5.32.0 or higher)

# Summary
## 1. Time series neural activity plot
## 2. Correlogram
## 3. Time series dimensionality reduction and clustering
## 4. Spectrogram
## 5. Time-delay embedding
## 6. Eigenvalue distribution of the history matrix of all cells
## 7. Recurrence plot

# How to reproduce this workflow

# Reference
1. [rikenbit/WormSCE](https://github.com/rikenbit/WormSCE)

# License
Copyright (c) 2020 Kentaro Yamamoto and RIKEN Bioinformatics Research Unit Released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

# Authors
- Kentaro Yamamoto
- Koki Tsuyuzaki
- Itoshi Nikaido