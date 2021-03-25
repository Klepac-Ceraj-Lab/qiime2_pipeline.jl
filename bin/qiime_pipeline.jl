#!/usr/bin julia

# This workflow was generated from LangilleLab workflow and Qiime Tutorials, Qiime version 2020.8 
# https://github.com/LangilleLab/microbiome_helper/wiki/Amplicon-SOP-v2-(qiime2-2020.8)

using ArgParse
using Logging
using LoggingExtras
using Dates
using QIIME2
using DataDeps

settings = ArgParseSettings()

@add_arg_table! settings begin
    "--verbose", "-v"
        help = "Set logging level to INFO"
        action = :store_true
    "--debug"
        help = "Set logging level to DEBUG. Overrides --verbose."
        action = :store_true
    "--quiet", "-q"
        help = "Supress all STDOUT logging other than ERROR. Overrides --verbose and --debug"
        action = :store_true
    "--log", "-l"
        help = "Store log to file, regardless of whether --quiet is set. That is, one can do `--verbose --log file.log --quiet`"
    "--threads"
        arg_type = Int
        help = "Set number of threads to use. Defaults to Threads.nthreads()"
        default = Threads.nthreads()

end

@add_arg_table! settings begin
    "basic"
        action = :command
end

# Project otions
@add_arg_table! settings["basic"] begin
    "--name", "-n"
        help = "Name of the project"
        default = "MyProject_$(today())"
    "--mapping_file", "-m"
        help = "Mapping sample names to metadata"
    "--classifier", "-c"
        help = "path to .qza file with classifier" 
        # default = datadep"classifier_138_99_full"
    "--input", "-i"
        help = "path to raw data input"
        default = "raw_data"
    "--output", "-o"
        help = "Directory for output. Defaults to current working directory"
        default = "./"
    "--force"
        help = "Re-run steps even if outputs directories already exist"
        action = :store_true
end

function basic(subargs)
    rawdir = abspath(subargs["input"])
    outdir = joinpath(abspath(subargs["output"]), subargs["name>`"])
    threads = get(subargs, "threads", Threads.nthreads())
    mapping_file = subargs["mapping_file"]
    
    @info subargs rawdir outdir threads mapping_file

    !isdir(outdir) && mkdir(outdir)

    fastqc_out = joinpath(outdir, "fastqc_out")
    if isdir(fastqc_out) && subargs["force"]
        @warn "Removing fastqc output dir: $fastqc_out"
        rm(fastqc_out, force = true)
    end

    if isdir(fastqc_out)
        @info "Fastqc directory exists, skipping step. Overwrite with --force"
    else
        isdir(subargs["input"]) || throw(ArgumentError("Input dir $(subargs["input"]) doesn't exist"))
        mkdir(fastqc_out)
        run(`fastqc $(readdir(rawdir, join=true)) -o $fastqc_out -t $threads`)
    end


    # Import Raw Sequence files as Qiime 2 artifact
    # This allows for standardization of QIIME 2 analyses and keeps track of all commands that were run to produce a file.
    # The extension for the artifact files is QZA.
    reads_qza = joinpath(outdir, "reads_qza")
    if isdir(reads_qza) && subargs["force"]
        @warn "Removing fastqc output dir: $reads_qza"
        rm(reads_qza, force = true)
    end

    if isdir(reads_qza)
        @info "Reads.qza directory exists, skipping import step. Overwrite with --force"
    else
        isdir(subargs["input"]) || throw(ArgumentError("Input dir $(subargs["input"]) doesn't exist"))
        mkdir(reads_qza)
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

    qiime_cmd("cutadapt", "trim-paired",
                i_demultiplexed_sequences = joinpath(reads_qza, "reads.qza"),
                p_cores                   = threads,
                p_front_f                 = "^GTGYCAGCMGCCGCGGTAA",
                p_front_r                 = "^CCGYCAATTYMTTTRAGTTT",
                o_trimmed_sequences       = joinpath(reads_qza, "reads_trimmed.qza")
                )

    dada2_out = joinpath(outdir, "dada2")
    if isdir(reads_qza) && subargs["force"]
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
    if isdir(reads_qza) && subargs["force"]
        @warn "Removing fastqc output dir: $taxa_out"
        rm(taxa_out, force = true)
    end

    if isdir(taxa_out)
        @info "Taxa directory exists, skipping classification step. Overwrite with --force"
    else
        classifier = subargs["classifier"]
        qiime_cmd("feature-classifier", "classify-sklearn";
                        i_reads = joinpath(dada2_out, "representative_sequences.qza"),
                        i_classifier = classifier,
                        p_n_jobs = threads,
                        output_dir = taxa_out
                    )

        qiime_cmd("tools", "export"; input_path = joinpath(taxa_out, "classification.qza"), output_path = taxa_out)
    end

    qiime_cmd("tools", "export",
                input_path = joinpath(taxa_out, "classification.qza"),
                output_path = taxa_dir
            )

    qiime_cmd("feature-table", "tabulate-seqs", 
                i_data = joinpath(dada2_out, "representative_sequences.qza"),
                o_visualization = joinpath(dada2_out, "representative_sequences.qzv")
            )

    # #!!!!change -p-min-frequency!!!
    # #filtering out rare ASVs (removed samples all samples that are <0.1% mean sample depth; mean sample depth =17,560 )
    qiime_cmd("feature-table", "filter-features",
                i_table = joinpath(dada2_out, "table.qza"),
                p_min_frequency = 20, # Should be 0.1 % of mean
                p_min_samples = 1,
                o_filtered_table = joinpath(dada2_out, "dada2_table_filt.qza")
            )

    qiime_cmd("taxa", "filter-table",
                i_table = joinpath(dada2_out, "dada2_table_filt.qza"),
                i_taxonomy = joinpath(taxa_out, "classification.qza"),
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
                p_metrics = "observed_features",
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

    qiime_cmd("taxa", "barplot",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                i_taxonomy = joinpath(taxa_out, "classification.qza"),
                m_metadata_file = mapping_file,
                o_visualization = joinpath(taxa_out, "taxa_barplot.qzv")
            )

    qiime_cmd("feature-table", "group",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                p_axis = "sample",
                p_mode = "sum",
                m_metadata_file = mapping_file,
                m_metadata_column = "FullTreatment",
                o_grouped_table = joinpath(dada2_out, "dada2_table_final_FullTreatment.qza")
            )

    diversity_dir = joinpath(outdir, "diversity")

    qiime_cmd("diversity", "core-metrics-phylogenetic",
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                i_phylogeny = joinpath(dada2_out, "rooted_tree.qza"),
                p_sampling_depth = 5000,
                m_metadata_file = mapping_file,
                output_dir = diversity_dir
            )

    qiime_cmd("diversity", "alpha-group-significance",
                i_alpha_diversity = joinpath(diversity_dir, "shannon_vector.qza"),
                m_metadata_file = mapping_file,
                o_visualization = joinpath(diversity_dir, "shannon_compare_groups.qzv")
            )

    qiime_cmd("diversity", "alpha-group-significance",
                i_alpha_diversity = joinpath(diversity_dir, "evenness_vector.qza"),
                m_metadata_file = mapping_file,
                o_visualization = joinpath(diversity_dir, "evenness_compare_groups.qzv")
            )

    qiime_cmd("diversity", "alpha-group-significance",
                i_alpha_diversity = joinpath(diversity_dir, "faith_pd_vector.qza"),
                m_metadata_file = mapping_file,
                o_visualization = joinpath(diversity_dir, "faith_pd_compare_groups.qzv")
            )

    qiime_cmd("composition", "add-pseudocount";
                i_table = joinpath(dada2_out, "dada2_table_final.qza"),
                p_pseudocount = 1,
                o_composition_table = joinpath(dada2_out, "dada2_table_final_pseudocount.qza")
            )

    ancom_dir = joinpath(outdir, "ancom")
    qiime_cmd("composition", "ancom",
                i_table = joinpath(dada2_out, "dada2_table_final_pseudocount.qza"),
                m_metadata_file = mapping_file,
                m_metadata_column = "FullTreatment",
                output_dir = ancom_dir
            )

    qiime_cmd("composition", "ancom",
                i_table = joinpath(dada2_out, "dada2_table_final_pseudocount.qza"),
                m_metadata_file = mapping_file,
                m_metadata_column = "Treatment",
                output_dir = ancom_dir*"2" # TODO - fix this
            )

    qiime_cmd("tools", "export",
                input_path = joinpath(dada2_out, "rep_seqs_final.qza"),
                output_path = "dada2_output_exported"
            )

    ## Need to figure out how to deal with duplicate subargs
    # qiime_cmd("longitudinal", "anova",
    #             m_metadata_file = joinpath(diversity_dir, "faith_pd_vector.qza"),
    #             m_metadata_file = mapping_file,
    #             p_formula = "'faith_pd ~ FullTreatment * Day'",
    #             o_visualization = joinpath(diversity_dir, "faiths_pd_anova.qzv"),
      
    # )

    qiime_cmd("diversity", "beta-group-significance",
                i_distance_matrix = joinpath(diversity_dir, "unweighted_unifrac_distance_matrix.qza"),
                m_metadata_file = mapping_file,
                m_metadata_column = "FullTreatment",
                o_visualization = joinpath(diversity_dir, "unweighted-unifrac-FullTreatment-significance.qzv"),

    )

    run(`sed -i -e '1 s/Feature/#Feature/' -e '1 s/Taxon/taxonomy/' taxa/taxonomy.tsv`)

    qiime_cmd("tools", "export",
                input_path = joinpath(dada2_out, "dada2_table_final.qza"),
                output_path = "dada2_output_exported"
    )


    biom_cmd("add-metadata";
            input_fp = "dada2_output_exported/feature-table.biom",
            output_fp = "dada2_output_exported/feature-table_w_tax.biom",
            observation_metadata_fp = joinpath(taxa_out, "taxonomy.tsv"),
            sc_separated = "taxonomy"
    )

    biom_cmd("convert";
            input_fp = "dada2_output_exported/feature-table_w_tax.biom",
            output_fp = "dada2_output_exported/feature-table_w_tax.txt",
            to_tsv = "",
            header_key = "taxonomy"
    )
end

function main()
    args = parse_args(ARGS, settings)

    logger = args["debug"] ? ConsoleLogger(stderr, Logging.Debug) :
             args["verbose"] ? MinLevelLogger(global_logger(), Logging.Info) :
             MinLevelLogger(global_logger(), Logging.Warn)

    flogger = !isnothing(args["log"]) ? MinLevelLogger(FileLogger(args["log"]), Logging.min_enabled_level(logger)) : nothing

    if args["quiet"]
        qlogger = ConsoleLogger(stderr, Logging.Error)
        isnothing(flogger) ? global_logger(qlogger) : global_logger(TeeLogger(qlogger, flogger))
    elseif isnothing(flogger)
        global_logger(logger)
    else
        global_logger(TeeLogger(logger, flogger))
    end

    torun = args["%COMMAND%"]

    Base.eval(@__MODULE__, Symbol(torun))(args[torun])
end

main()

main()
