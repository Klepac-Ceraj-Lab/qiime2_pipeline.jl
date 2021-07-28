module QIIME2

using DataDeps
using CSV
using DataFrames
using Statistics
# using Conda

export 
    qiime_cmd,
    biom_cmd,
    depth_filter

function __init__()
    register(DataDep(
        "classifier_138_99_full",
        """
        Silva 138 99% OTUs full-length sequences (MD5: fddefff8bfa2bbfa08b9cad36bcdf709)
        """,
        "https://data.qiime2.org/2021.2/common/silva-138-99-nb-classifier.qza",
        "def48c9f9c8c3444f42b13dbeaf5f6376efff3e8e81994788dc3493fe02aaedc" # md5sum fddefff8bfa2bbfa08b9cad36bcdf709
        )
    )
end

function qiime_cmd(command, subcommand; kwargs...)
    args = Iterators.flatten((replace(string("--", k), "_"=>"-"), v) for (k,v) in pairs(kwargs))
    cmd = `qiime $command $subcommand $args`
    @info "Running command: $cmd"
    return run(cmd)
end

function biom_cmd(command; kwargs...)
    c = ["biom", command]
    append!(c, Iterators.flatten((replace(string("--", k), "_"=>"-"), v) for (k,v) in pairs(kwargs)))
    deleteat!(c, findall(==(""), c))
    @info "Running command: $(Cmd(c))"
    return run(Cmd(c))
end


function depth_filter(dada_dir, perc=0.1)
    biom = joinpath(dada_dir, "feature-table.biom")
    features = joinpath(dada_dir, "feature-table.tsv")
    
    if !isfile(features)
        @info "Couldn't find feature table, trying to export from qza form"
        qiime_cmd("tools", "export",
                   input_path = joinpath(dada_dir, "table.qza"),
                   output_path = dada_dir
        )
        biom_cmd("convert", input_fp=biom, output_fp=features, to_tsv="")
    end

    df = CSV.read(features, DataFrame, header=2, delim='\t')
    return mean(sum.(eachcol(df[!, 2:end]))) * (perc / 100)
end

end