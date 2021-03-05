#!/usr/bin julia

# This workflow was generated from LangilleLab workflow and Qiime Tutorials, Qiime version 2020.8 
# https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2020.8)

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using ArgParse
using Logging
using LoggingExtras
using Dates
using Threads
using QIIME2

s = ArgParseSettings()

# Project otions
@add_arg_table! s begin
    "project_name"
        help = "Name of the project"
        default = "MyProject_$(today())"
    "mapping_file"
        help = "Mapping sample names to metadata"
    "classifier" # TODO: get from datadeps
        help = "path to .qza file with classifier" 
        default = datadep"classifier_138_99_full/silva-138-99-nb-classifier.qza"

    "--input", "-i"
        help = "path to raw data input"
        default = "raw_data"
    "--output", "-o"
        help = "Directory for output. Defaults to current working directory"
        default = "./"
    "--threads"
        type = Int
        help = "Set number of threads to use. Defaults to Threads.nthreads()"
    "--force"
        help = "Re-run steps even if outputs directories already exist"
        action = :store_true
end

# Logging options
@add_arg_table! s begin
    "--verbose", "-v"
        help = "Set logging level to INFO"
        action = :store_true
    "--debug"
        help = "Set logging level to DEBUG"
        action = :store_true
    "--quiet", "-q"
        help = "Supress STDOUT logging"
        action = :store_true
end


function main()
    args = parse_args(ARGS, s)
    rawdir = abspath(args["input"])
    outdir = abspath(args["output"])
    threads = get(args, "threads", Threads.nthreads())

    fastqc_out = joinpath(outdir, "fastqc_out")
    if isdir(fastqc_out) && args["force"]
        @warn "Removing fastqc output dir: $fastqc_out"
        rm(fastqc_out, force = true)
    end

    if isdir(fastqc_out)
        @info "Fastqc directory exists, skipping step. Overwrite with --force"
    else
        isdir(args["input"]) || throw(ArgumentError("Input dir $(args["input"]) doesn't exist"))
        run(`fastqc $rawdir -o $fastqc_out -t $threads`)
    end


    # Import Raw Sequence files as Qiime 2 artifact
    # This allows for standardization of QIIME 2 analyses and keeps track of all commands that were run to produce a file.
    # The extension for the artifact files is QZA.
    reads_qza = joinpath(outdir, "reads_qza")
    if isdir(reads_qza) && args["force"]
        @warn "Removing fastqc output dir: $reads_qza"
        rm(reads_qza, force = true)
    end

    if isdir(reads_qza)
        @info "Reads.qza directory exists, skipping import step. Overwrite with --force"
    else
        isdir(args["input"]) || throw(ArgumentError("Input dir $(args["input"]) doesn't exist"))
        qiime_cmd("tools", "import";
                  type         = "SampleData[PairedEndSequencesWithQuality]",
                  input_path   = rawdir,
                  output_path  = joinpath(reads_qza, "reads.qza"),
                  input_format = "CasavaOneEightSingleLanePerSampleDirFmt" # TODO: Make input format configurable
        ) 
    end

    # Trim primers with cutadapt
    # Trim primers with cutadapt - remove all primers and take out all sequences that don't begin with primer sequence  
    # 16S V4-V5 region bacteria + archaea  
    # 515F = GTGYCAGCMGCCGCGGTAA  
    # 926R = CCGYCAATTYMTTTRAGTTT  

    # For V6-V8 region:  
    #    --p-front-f TYAATYGGANTCAACRCC
    #    --p-front-r CRGTGWGTRCAAGGRGCA

    if isdir(reads_qza)
        @info "Reads.qza directory exists, skipping cutadapt step. Overwrite with --force"
    else
        qiime_cmd("cutadapt", "trim-paired",
                    i_demultiplexed_sequences = joinpath(reads_qza, "reads.qza"),
                    p_cores                   = threads,
                    p_front_f                 = "^GTGYCAGCMGCCGCGGTAA",
                    p_front_r                 = "^CCGYCAATTYMTTTRAGTTT",
                    o_trimmed_sequences       = joinpath(reads_qza, "reads_trimmed.qza")
                    )
    end

    dada2_out = joinpath(outdir, "dada2")
    if isdir(reads_qza) && args["force"]
        @warn "Removing fastqc output dir: $dada2_out"
        rm(dada2_out, force = true)
    end

    if isdir(dada2_out)
        @info "Dada2 directory exists, skipping step. Overwrite with --force"
    else
        # TODO: Make options configurable
        qiime_cmd("dada2", "denoise-paired";
                        i_demultiplexed_seqs = joinpath(reads_qza, "reads_trimmed.qza"),
                        p_trunc_len_f        = 270,
                        p_trunc_len_r        = 210,
                        p_max_ee_f           = 2,
                        p_max_ee_r           = 3,
                        p_n_threads          = threads,
                        output_dir           = dada2_out
                    ) 
    end

    taxa_out = joinpath(outdir, "taxa")
    if isdir(reads_qza) && args["force"]
        @warn "Removing fastqc output dir: $taxa_out"
        rm(taxa_out, force = true)
    end

    if isdir(taxa_out)
        @info "Taxa directory exists, skipping classification step. Overwrite with --force"
    else
        classifier = args["classifier"]
        qiime_cmd("feature-classifier", "classify-sklearn";
                        i_reads = joinpath(dada2_out, "representative_sequences.qza"),
                        i_classifier = classifier,
                        p_n_jobs = threads,
                        output_dir = taxa_out
                    )

        qiime_cmd("tools", "export"; input_path = joinpath(taxa_out, "classification.qza"), output_path = taxa_out)
    end

    qiime_cmd("tools", "export",
                input_path = joinpath("taxa_out", "classification.qza"),
                output_path = "taxa"
            )

    qiime_cmd("feature-table", "tabulate-seqs", 
                i_data = joinpath(dada2_out, "representative_sequences.qza"),
                o_visualization = joinpath(dada2_out, "representative_sequences.qzv")
            )

    # #!!!!change -p-min-frequency!!!
    # #filtering out rare ASVs (removed samples all samples that are <0.1% mean sample depth; mean sample depth =17,560 )
    qiime_cmd("feature-table", "filter-features",
                i_table = "dada2_output/table.qza",
                p_min_frequency = 20,
                p_min_samples = 1,
                o_filtered_table = joinpath(dada2_out, "dada2_table_filt.qza")
            )

    qiime_cmd("taxa", "filter-table",
                i_table = joinpath(dada2_out, "dada2_table_filt.qza"),
                i_taxonomy = joinpath("taxa_out", "classification.qza"),
                p_include = "p__",
                p_exclude = "mitochondria,chloroplast",
                o_filtered_table = joinpath(dada2_out, "dada2_table_filt_contam.qza")   
            )


    cp(joinpath(dada2_out, "dada2_table_filt_contam.qza"), joinpath(dada2_out, "dada2_table_final.qza"))

    #p-max-depth needs to be changed
    qiime_cmd("diversity", "alpha-rarefaction",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                p_max_depth = 36000,
                p_steps = 20,
                p_metrics = "'observed_features'",
                o_visualization = joinpath(dada2_out, "rarefaction_curves_test.qzv")
            )

    qiime_cmd("feature-table", "filter-samples",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                p_min_frequency = 5000,
                o_filtered_table = joinpath(dada2_out, "dada2_table_final.qza")
            )

    qiime_cmd("feature-table", "filter-seqs",
                i_data = joinpath(dada2_out, "representative_sequences.qza"),
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                o_filtered_data = joinpath(dada2_out, "rep_seqs_final.qza")
            )

    qiime_cmd("feature-table", "summarize",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                o_visualization = joinpath(dada2_out, "dada2_table_final_summary.qzv")
            )

    qiime_cmd("phylogeny", "align-to-tree-mafft-fasttree",
                i_sequences = joinpath(dada2_out, "rep_seqs_final.qza"),
                o_alignment = joinpath(dada2_out, "aligned_rep_seqs.qza"),
                o_masked_alignment = joinpath(dada2_out, "masked_aligned_rep_seqs.qza"),
                o_tree = joinpath(dada2_out, "unrooted_tree.qza"),
                o_rooted_tree = joinpath(dada2_out, "rooted_tree.qza")
            )

    qiime_cmd("diversity", "alpha-rarefaction",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                p_max_depth = 36000,
                p_steps = 20,
                i_phylogeny = joinpath(dada2_out, "rooted_tree.qza"),
                m_metadata_file = mapping_file,
                o_visualization = "rarefaction_curves.qzv"
            )

    qiime_cmd("diversity", "alpha-rarefaction",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                p_max_depth = 36000,
                p_steps = 20,
                i_phylogeny = joinpath(dada2_out, "rooted_tree.qza"),
                o_visualization = "rarefaction_curves_eachsample.qzv"
            )

    qiime_cmd("taxa", "barplot",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                i_taxonomy = joinpath("taxa_out", "classification.qza"),
                m_metadata_file = mapping_file,
                o_visualization = joinpath("taxa_out", "taxa_barplot.qzv")
            )

    qiime_cmd("feature-table", "group",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                p_axis = "sample",
                p_mode = "sum",
                m_metadata_file = mapping_file,
                m_metadata_column = "FullTreatment",
                o_grouped_table = joinpath(dada2_out, "dada2_table_final_FullTreatment.qza")
            )

    qiime_cmd("diversity", "core-metrics-phylogenetic",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                i_phylogeny = joinpath(dada2_out, "rooted_tree.qza"),
                p_sampling_depth = 5000,
                m_metadata_file = mapping_file,
                output_dir = "diversity"
            )

    qiime_cmd("diversity", "alpha-group-significance",
                i_alpha_diversity = "diversity/shannon_vector.qza",
                m_metadata_file = mapping_file,
                o_visualization = "diversity/shannon_compare_groups.qzv"
            )

    qiime_cmd("diversity", "alpha-group-significance",
                i_alpha_diversity = "diversity/evenness_vector.qza",
                m_metadata_file = mapping_file,
                o_visualization = "diversity/evenness_compare_groups.qzv"
            )

    qiime_cmd("diversity", "alpha-group-significance",
                i_alpha_diversity = "diversity/faith_pd_vector.qza",
                m_metadata_file = mapping_file,
                o_visualization = "diversity/faith_pd_compare_groups.qzv"
            )

    qiime_cmd("composition", "add-pseudocount",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                p_pseudocount = 1,
                o_composition_table = joinpath(dada2_out, "dada2_table_final_pseudocount.qza")
            )

    qiime_cmd("composition", "ancom",
                i_table = joinpath(dada2_out, "dada2_table_final_pseudocount.qza"),
                m_metadata_file = mapping_file,
                m_metadata_column = "FullTreatment",
                output_dir = "ancom_output"
            )

    qiime_cmd("composition", "ancom",
                i_table = joinpath(dada2_out, "dada2_table_final_pseudocount.qza"),
                m_metadata_file = mapping_file,
                m_metadata_column = "Treatment",
                output_dir = "ancom2_output    "
            )

    qiime_cmd("tools", "export",
                input_path = joinpath(dada2_out, "rep_seqs_final.qza"),
                output_path = "dada2_output_exported"
            )

    qiime_cmd("longitudinal", "anova",
                m_metadata_file = "diversity/faith_pd_vector.qza",
                m_metadata_file = mapping_file,
                p_formula = "'faith_pd ~ FullTreatment * Day'",
                o_visualization = "diversity/faiths_pd_anova.qzv",
      
    )

    qiime_cmd("diversity", "beta-group-significance",
                i_distance_matrix = "diversity/unweighted_unifrac_distance_matrix.qza",
                m_metadata_file = mapping_file,
                m_metadata_column = "FullTreatment",
                o_visualization = "diversity/unweighted-unifrac-FullTreatment-significance.qzv",

    )

    qiime_cmd("diversity", "beta-group-significance",
                i_distance_matrix = "diversity/weighted_unifrac_distance_matrix.qza",
                m_metadata_file = mapping_file,
                m_metadata_column = "FullTreatment",
                o_visualization = "diversity/weighted-unifrac-FullTreatment-significance.qzv",

    )

    qiime_cmd("diversity", "beta-group-significance",
                i_distance_matrix = "diversity/bray_curtis_distance_matrix.qza",
                m_metadata_file = mapping_file,
                m_metadata_column = "FullTreatment",
                o_visualization = "diversity/Bray-Curtis-FullTreatment-significance.qzv",

    )

    run(`sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv`)

    qiime_cmd("tools", "export",
                input_path = joinpath(dada2_out, "dada2_table_final.qza"),
                output_path = "dada2_output_exported"
    )


    biom_cmd("add-metadata";
            input = "dada2_output_exported/feature-table.biom",
            output = "dada2_output_exported/feature-table_w_tax.biom",
            observation_metadata_fp = "taxa/taxonomy.tsv",
            sc_separated = "taxonomy"
    )

    biom_cmd("convert";
            input = "dada2_output_exported/feature-table_w_tax.biom",
            output = "dada2_output_exported/feature-table_w_tax.txt",
            to_tsv = "header_key = \"taxonomy\""
    )


end
main()



# ```python
# !



# ```

#     sed: -e: No such file or directory
#     [32mExported dada2_output/dada2_table_final.qza as BIOMV210DirFmt to directory dada2_output_exported[0m
#     Traceback (most recent call last):
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/bin/biom", line 11, in <module>
#         sys.exit(cli())
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 829, in __call__
#         return self.main(*args, **kwargs)
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 782, in main
#         rv = self.invoke(ctx)
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 1259, in invoke
#         return _process_result(sub_ctx.command.invoke(sub_ctx))
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 1066, in invoke
#         return ctx.invoke(self.callback, **ctx.params)
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/click/core.py", line 610, in invoke
#         return callback(*args, **kwargs)
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/biom/cli/metadata_adder.py", line 107, in add_metadata
#         float_fields, sample_header, observation_header)
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/biom/cli/metadata_adder.py", line 174, in _add_metadata
#         header=observation_header)
#       File "/Users/vklepacc/miniconda3/envs/qiime2-2020.8/lib/python3.6/site-packages/biom/parse.py", line 538, in from_file
#         raise BiomParseException("No header line was found in mapping "
#     biom.exception.BiomParseException: No header line was found in mapping file.
#     Usage: biom convert [OPTIONS]
#     Try 'biom convert -h' for help.
    
#     Error: Invalid value for '-i' / '--input-fp': File 'dada2_output_exported/feature-table_w_tax.biom' does not exist.
