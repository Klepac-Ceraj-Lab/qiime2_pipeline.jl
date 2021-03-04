# Nisreen's Bokashi Project -- 16S rRNA gene data
This workflow was generated from LangilleLab workflow and Qiime Tutorials, Qiime version 2020.8 

https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2020.8)

In terminal:
conda source activate qiime2-2020.8
jupyter notebook


```python
from os import chdir, mkdir
from os.path import join
#the following are only available in the current development branch of IPython
from IPython.display import FileLinks, FileLink

%pylab inline
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Set some Pandas options
pd.set_option('display.notebook_repr_html', False)
pd.set_option('display.max_columns', 20)
pd.set_option('display.max_rows', 25)
#populating the interactive namespace for numpy and matplotlib - not needed now...
```

    Populating the interactive namespace from numpy and matplotlib



```python
project_name = "Nisreen_bokashi"
mapping_file = "./Nisreen_16Smapping.tsv"
classifier = "/Users/vklepacc/classifiers/silva-138-99-nb-classifier.qza"
NCORES = 2
#to avoid changing mapping file / classifier throughout the code, just refer to it by typing $mapping_file
```

## Inspect read quality
Need FastQC and MultiQC (separately downloaded as these are not packaged in qiime2)


```python
    !mkdir fastqc_out
    !fastqc -t 4 raw_data/* -o fastqc_out
```


      File "<ipython-input-4-c5a0fe8269cf>", line 2
        fastqc -t 4 raw_data/* -o fastqc_out
                  ^
    SyntaxError: invalid syntax



## Import Raw Sequence files as Qiime 2 artifact
This allows for standardization of QIIME 2 analyses and keeps track of all commands that were run to produce a file. The extension for the artifact files is QZA.


```python
!mkdir reads_qza

!qiime tools import --type SampleData[PairedEndSequencesWithQuality] \
                   --input-path raw_data \
                   --output-path reads_qza/reads.qza \
                   --input-format CasavaOneEightSingleLanePerSampleDirFmt
```

    [32mImported raw_data as CasavaOneEightSingleLanePerSampleDirFmt to reads_qza/reads.qza[0m


## Trim primers with cutadapt
Trim primers with cutadapt - remove all primers and take out all sequences that don't begin with primer sequence  
16S V4-V5 region bacteria + archaea  
515F = GTGYCAGCMGCCGCGGTAA  
926R = CCGYCAATTYMTTTRAGTTT  

For V6-V8 region:  
   --p-front-f TYAATYGGANTCAACRCC \
   --p-front-r CRGTGWGTRCAAGGRGCA \


```python
!qiime cutadapt trim-paired --i-demultiplexed-sequences reads_qza/reads.qza \
                           --p-cores $NCORES \
                           --p-front-f ^GTGYCAGCMGCCGCGGTAA \
                           --p-front-r ^CCGYCAATTYMTTTRAGTTT \
                           --o-trimmed-sequences reads_qza/reads_trimmed.qza
```

    [32mSaved SampleData[PairedEndSequencesWithQuality] to: reads_qza/reads_trimmed.qza[0m


## Summarize trimmed FASTQs
qiime demux summarize command gives back a report of the number of reads per sample and quality distribution across the reads.


```python
!qiime demux summarize \
   --i-data reads_qza/reads_trimmed.qza \
   --o-visualization reads_qza/reads_trimmed_summary.qzv
```

    [32mSaved Visualization to: reads_qza/reads_trimmed_summary.qzv[0m


## Denoising reads into ASVs
DADA2 - https://benjjneb.github.io/dada2/tutorial.html


```python
#!qiime dada2 denoise-paired --help
```


```python
!qiime dada2 denoise-paired --i-demultiplexed-seqs reads_qza/reads_trimmed.qza \
                          --p-trunc-len-f 270 \
                          --p-trunc-len-r 210 \
                          --p-max-ee-f 2 \
                          --p-max-ee-r 3 \
                          --p-n-threads 1 \ 
                          --output-dir dada2_output --verbose
#p-n-threads should prob. be higher on a server 
```

    Running external command line application(s). This may print messages to stdout and/or stderr.
    The command(s) being run are below. These commands cannot be manually re-run as they will depend on temporary files that no longer exist.
    
    Command: run_dada_paired.R /var/folders/k0/my9vv5l15_dfcwf1xl3dksg00000gp/T/tmpenxz4c2i/forward /var/folders/k0/my9vv5l15_dfcwf1xl3dksg00000gp/T/tmpenxz4c2i/reverse /var/folders/k0/my9vv5l15_dfcwf1xl3dksg00000gp/T/tmpenxz4c2i/output.tsv.biom /var/folders/k0/my9vv5l15_dfcwf1xl3dksg00000gp/T/tmpenxz4c2i/track.tsv /var/folders/k0/my9vv5l15_dfcwf1xl3dksg00000gp/T/tmpenxz4c2i/filt_f /var/folders/k0/my9vv5l15_dfcwf1xl3dksg00000gp/T/tmpenxz4c2i/filt_r 270 210 0 0 2.0 3.0 2 independent consensus 1.0 1 1000000
    
    R version 3.5.1 (2018-07-02) 
    Loading required package: Rcpp
    DADA2: 1.10.0 / Rcpp: 1.0.4.6 / RcppParallel: 5.0.0 
    1) Filtering ......................................................................................
    2) Learning Error Rates
    277396920 total bases in 1027396 reads from 37 samples will be used for learning the error rates.
    215753160 total bases in 1027396 reads from 37 samples will be used for learning the error rates.
    3) Denoise samples ......................................................................................
    ......................................................................................
    4) Remove chimeras (method = consensus)
    6) Write output
    [32mSaved FeatureTable[Frequency] to: dada2_output/table.qza[0m
    [32mSaved FeatureData[Sequence] to: dada2_output/representative_sequences.qza[0m
    [32mSaved SampleData[DADA2Stats] to: dada2_output/denoising_stats.qza[0m



```python
#for checking how many reads are retained
!qiime tools export --input-path dada2_output/denoising_stats.qza --output-path dada2_output
```

    [32mExported dada2_output/denoising_stats.qza as DADA2StatsDirFmt to directory dada2_output[0m



```python
#summarizing DADA2 output
!qiime feature-table summarize --i-table dada2_output/table.qza --o-visualization dada2_output/table_summary.qzv
```

    [32mSaved Visualization to: dada2_output/table_summary.qzv[0m


# Assign taxonomy to ASVs
You can assign taxonomy to your ASVs using a Naive-Bayes approach implemented in the scikit learn Python library and the SILVA database. 
  
### Build or acquire taxonomic classifier
The full-length 16S/18S classifier was downloaded from the QIIME 2 website (silva-138-99-nb-classifier.qza for the latest classifier). 
 
### Run taxonomic classification
You can run the taxonomic classification with this command, which is one of the longest running and most memory-intensive command of the SOP.


```python
#!qiime feature-classifier classify-sklearn --help
```


```python
!qiime feature-classifier classify-sklearn --i-reads dada2_output/representative_sequences.qza\
                                          --i-classifier $classifier \
                                          --p-n-jobs $NCORES \
                                          --output-dir taxa \
--verbose

```

    [32mSaved FeatureData[Taxonomy] to: taxa/classification.qza[0m


As with all QZA files, you can export the output file to take a look at the classifications and confidence scores:


```python
!qiime tools export --input-path taxa/classification.qza --output-path taxa
```

    [32mExported taxa/classification.qza as TSVTaxonomyDirectoryFormat to directory taxa[0m


## Assess subset of taxonomic scores by blast
The performance of the taxonomic classification is difficult to assess without a gold-standard reference, but nonetheless one basic sanity check is to compare the taxonomic assignments with the top BLASTn hits for certain ASVs.



```python
!qiime feature-table tabulate-seqs --i-data dada2_output/representative_sequences.qza \
                                   --o-visualization dada2_output/representative_sequences.qzv
```

    [32mSaved Visualization to: dada2_output/representative_sequences.qzv[0m


Fiter the representative_sequences.qza table  
- 


```python
#!!!!change -p-min-frequency!!!
#filtering out rare ASVs (removed samples all samples that are <0.1% mean sample depth; mean sample depth =17,560 )
!qiime feature-table filter-features \
   --i-table dada2_output/table.qza \
   --p-min-frequency 20 \
   --p-min-samples 1 \
   --o-filtered-table dada2_output/dada2_table_filt.qza
```

    [32mSaved FeatureTable[Frequency] to: dada2_output/dada2_table_filt.qza[0m



```python
#filtering out contaminant and unclassified ASVs
!qiime taxa filter-table \
   --i-table dada2_output/dada2_table_filt.qza \
   --i-taxonomy taxa/classification.qza \
   --p-include p__ \
   --p-exclude mitochondria,chloroplast \
   --o-filtered-table dada2_output/dada2_table_filt_contam.qza
```

    [32mSaved FeatureTable[Frequency] to: dada2_output/dada2_table_filt_contam.qza[0m



```python
!qiime feature-table summarize \
   --i-table dada2_output/dada2_table_filt_contam.qza \
   --o-visualization dada2_output/dada2_table_filt_contam_summary.qzv
```

    [32mSaved Visualization to: dada2_output/dada2_table_filt_contam_summary.qzv[0m



```python
cp dada2_output/dada2_table_filt_contam.qza dada2_output/dada2_table_final.qza
```


```python
#p-max-depth needs to be changed
!qiime diversity alpha-rarefaction \
   --i-table dada2_output/dada2_table_final.qza \
   --p-max-depth 36000 \
   --p-steps 20 \
   --p-metrics 'observed_features' \
   --o-visualization dada2_output/rarefaction_curves_test.qzv
```

    [32mSaved Visualization to: dada2_output/rarefaction_curves_test.qzv[0m



```python
!qiime feature-table filter-samples \
   --i-table dada2_output/dada2_table_final.qza \
   --p-min-frequency 5000 \
   --o-filtered-table dada2_output/dada2_table_final.qza
#setting a cutoff at 5000, BI4-D3 did not make it to the list. CF did not make it either
```

    [32mSaved FeatureTable[Frequency] to: dada2_output/dada2_table_final.qza[0m


## filter tables and sequence data


```python
!qiime feature-table filter-seqs \
   --i-data dada2_output/representative_sequences.qza \
   --i-table dada2_output/dada2_table_final.qza \
   --o-filtered-data dada2_output/rep_seqs_final.qza

!qiime feature-table summarize \
   --i-table dada2_output/dada2_table_final.qza \
   --o-visualization dada2_output/dada2_table_final_summary.qzv
```

    [32mSaved FeatureData[Sequence] to: dada2_output/rep_seqs_final.qza[0m
    [32mSaved Visualization to: dada2_output/dada2_table_final_summary.qzv[0m


## build tree
there's a newer way to add sequences to a tree that I haven't yet used and it uses the sepp command. looks better bc the short reads are already added to the existing tree.


```python
!qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences dada2_output/rep_seqs_final.qza \
  --o-alignment dada2_output/aligned-rep-seqs.qza \
  --o-masked-alignment dada2_output/masked-aligned-rep-seqs.qza \
  --o-tree dada2_output/unrooted-tree.qza \
  --o-rooted-tree dada2_output/rooted-tree.qza
```

    [32mSaved FeatureData[AlignedSequence] to: dada2_output/aligned-rep-seqs.qza[0m
    [32mSaved FeatureData[AlignedSequence] to: dada2_output/masked-aligned-rep-seqs.qza[0m
    [32mSaved Phylogeny[Unrooted] to: dada2_output/unrooted-tree.qza[0m
    [32mSaved Phylogeny[Rooted] to: dada2_output/rooted-tree.qza[0m



```python
#generate rarefaction curves
!qiime diversity alpha-rarefaction \
   --i-table dada2_output/dada2_table_final.qza \
   --p-max-depth 36000 \
   --p-steps 20 \
   --i-phylogeny dada2_output/rooted-tree.qza \
   --m-metadata-file $mapping_file \
   --o-visualization rarefaction_curves.qzv
```

    [32mSaved Visualization to: rarefaction_curves.qzv[0m



```python
!qiime diversity alpha-rarefaction \
   --i-table dada2_output/dada2_table_final.qza \
   --p-max-depth 36000 \
   --p-steps 20 \
   --i-phylogeny dada2_output/rooted-tree.qza \
   --o-visualization rarefaction_curves_eachsample.qzv
```


```python
#generate stacked barchart of relative taxon abundances
!qiime taxa barplot \
   --i-table dada2_output/dada2_table_final.qza \
   --i-taxonomy taxa/classification.qza \
   --m-metadata-file $mapping_file \
   --o-visualization taxa/taxa_barplot.qzv
```


```python
!qiime feature-table group \
   --i-table dada2_output/dada2_table_final.qza \
   --p-axis sample \
   --p-mode sum \
   --m-metadata-file $mapping_file \
   --m-metadata-column FullTreatment \
   --o-grouped-table dada2_output/dada2_table_final_FullTreatment.qza
```

Calculating diversity metrics and generating ordination plots
Common alpha and beta-diversity metrics can be calculated with a single command in QIIME2. In addition, ordination plots (such as PCoA plots for weighted UniFrac distances) will be generated automatically as well. This command will also rarefy all samples to the sample sequencing depth before calculating these metrics (X is a placeholder for the lowest reasonable sample depth; samples with depth below this cut-off will be excluded).


```python
!qiime diversity core-metrics-phylogenetic \
   --i-table dada2_output/dada2_table_final.qza \
   --i-phylogeny dada2_output/rooted-tree.qza \
   --p-sampling-depth 5000 \
   --m-metadata-file $mapping_file  \
   --output-dir diversity
```


```python
!qiime diversity alpha-group-significance \
   --i-alpha-diversity diversity/shannon_vector.qza \
   --m-metadata-file $mapping_file \
   --o-visualization diversity/shannon_compare_groups.qzv
```


```python
!qiime diversity alpha-group-significance \
   --i-alpha-diversity diversity/evenness_vector.qza \
   --m-metadata-file $mapping_file \
   --o-visualization diversity/evenness_compare_groups.qzv

!qiime diversity alpha-group-significance \
   --i-alpha-diversity diversity/faith_pd_vector.qza \
   --m-metadata-file $mapping_file \
   --o-visualization diversity/faith_pd_compare_groups.qzv
```


```python
#add pseudocount... cannot have non-zero abundances
!qiime composition add-pseudocount \
   --i-table dada2_output/dada2_table_final.qza \
   --p-pseudocount 1 \
   --o-composition-table dada2_output/dada2_table_final_pseudocount.qza
```


```python
!qiime composition ancom \
   --i-table dada2_output/dada2_table_final_pseudocount.qza \
   --m-metadata-file $mapping_file \
   --m-metadata-column FullTreatment \
   --output-dir ancom_output
```


```python
!qiime composition ancom \
   --i-table dada2_output/dada2_table_final_pseudocount.qza \
   --m-metadata-file $mapping_file \
   --m-metadata-column Treatment \
   --output-dir ancom2_output
```


```python
!qiime tools export \
   --input-path dada2_output/rep_seqs_final.qza \
   --output-path dada2_output_exported
```

    [32mExported dada2_output/rep_seqs_final.qza as DNASequencesDirectoryFormat to directory dada2_output_exported[0m



```python
!qiime longitudinal anova \
  --m-metadata-file diversity/faith_pd_vector.qza \
  --m-metadata-file $mapping_file \
  --p-formula 'faith_pd ~ FullTreatment * Day' \
  --o-visualization diversity/faiths_pd_anova.qzv
```

    [32mSaved Visualization to: diversity/faiths_pd_anova.qzv[0m



```python
!qiime diversity beta-group-significance \
  --i-distance-matrix diversity/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $mapping_file \
  --m-metadata-column FullTreatment \
  --o-visualization diversity/unweighted-unifrac-FullTreatment-significance.qzv

!qiime diversity beta-group-significance \
  --i-distance-matrix diversity/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file $mapping_file \
  --m-metadata-column FullTreatment \
  --o-visualization diversity/weighted-unifrac-FullTreatment-significance.qzv

!qiime diversity beta-group-significance \
  --i-distance-matrix diversity/bray_curtis_distance_matrix.qza \
  --m-metadata-file $mapping_file \
  --m-metadata-column FullTreatment \
  --o-visualization diversity/Bray-Curtis-FullTreatment-significance.qzv
```

    [32mSaved Visualization to: diversity/unweighted-unifrac-FullTreatment-significance.qzv[0m
    [32mSaved Visualization to: diversity/weighted-unifrac-FullTreatment-significance.qzv[0m
    [32mSaved Visualization to: diversity/Bray-Curtis-FullTreatment-significance.qzv[0m



```python
!qiime longitudinal volatility \
  --m-metadata-file $mapping_file \
  --m-metadata-file diversity/unweighted_unifrac_pcoa_results.qza \
  --p-state-column Time \
  --p-individual-id-column Sample \
  --p-default-group-column 'FullTreatment' \
  --p-default-metric 'Axis 2' \
  --o-visualization ./pc_vol.qzv
```

    [32mSaved Visualization to: ./pc_vol.qzv[0m



```python
!sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv

!qiime tools export \
   --input-path dada2_output/dada2_table_final.qza \
   --output-path dada2_output_exported

!biom add-metadata \
   -i dada2_output_exported/feature-table.biom \
   -o dada2_output_exported/feature-table_w_tax.biom \
   --observation-metadata-fp taxa/taxonomy.tsv \
   --sc-separated taxonomy

!biom convert \
   -i dada2_output_exported/feature-table_w_tax.biom \
   -o dada2_output_exported/feature-table_w_tax.txt \
   --to-tsv \
   --header-key taxonomy
```

    sed: -e: No such file or directory
    [32mExported dada2_output/dada2_table_final.qza as BIOMV210DirFmt to directory dada2_output_exported[0m
    Traceback (most recent call last):
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/bin/biom", line 11, in <module>
        sys.exit(cli())
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 829, in __call__
        return self.main(*args, **kwargs)
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 782, in main
        rv = self.invoke(ctx)
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 1259, in invoke
        return _process_result(sub_ctx.command.invoke(sub_ctx))
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 1066, in invoke
        return ctx.invoke(self.callback, **ctx.params)
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 610, in invoke
        return callback(*args, **kwargs)
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/biom/cli/metadata_adder.py", line 107, in add_metadata
        float_fields, sample_header, observation_header)
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/biom/cli/metadata_adder.py", line 174, in _add_metadata
        header=observation_header)
      File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/biom/parse.py", line 538, in from_file
        raise BiomParseException("No header line was found in mapping "
    biom.exception.BiomParseException: No header line was found in mapping file.
    Usage: biom convert [OPTIONS]
    Try 'biom convert -h' for help.
    
    Error: Invalid value for '-i' / '--input-fp': File 'dada2_output_exported/feature-table_w_tax.biom' does not exist.



```python

```
