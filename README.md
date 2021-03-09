

## Conda environment

```sh
$ conda create -n qiime python=3.6 -y # requied by qiime2
$ conda activate qiime
$ conda install -yc qiime2 qiime2 q2cli
$ conda install -y fastqc
$ conda install -y cutadapt
$ conda install -y biom
$ conda install -c qiime2 q2-alignment q2-alignment q2-cutadapt q2-dada2 q2-taxa q2-feature-table q2-phylogeny q2-composition q2-feature-classifier q2-diversity
```