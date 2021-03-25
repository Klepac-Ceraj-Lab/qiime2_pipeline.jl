module QIIME2

using DataDeps
# using Conda

export 
    qiime_cmd,
    biom_cmd

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

end