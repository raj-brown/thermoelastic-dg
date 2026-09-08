#!/usr/bin/env julia

const REPO_ROOT = @__DIR__

const CASES = Dict(
    "case_1" => ("Carcione Figure 6 heat-source experiment", "case_1_carcione.jl"),
    "case_2" => ("Carcione Figure 7 comparison experiment", "case_2_carcione.jl"),
)

function print_usage()
    println("Usage: julia --project=. main.jl [all|case_1|case_2|--list]")
    println()
    println("Cases:")
    for name in sort!(collect(keys(CASES)))
        description, script = CASES[name]
        println("  $(rpad(name, 8)) $(description) ($(script))")
    end
    println("  all      Run every case in order (default)")
end

function run_case(name::String)
    description, script = CASES[name]
    path = joinpath(REPO_ROOT, script)
    isfile(path) || error("Case script not found: $path")

    println()
    println("=" ^ 78)
    println("Running $name: $description")
    println("Script: $path")
    println("=" ^ 78)

    run(`$(Base.julia_cmd()) --project=$(REPO_ROOT) $path`)
end

function main(args=ARGS)
    if any(arg -> arg in ("-h", "--help", "--list"), args)
        print_usage()
        return nothing
    end

    isempty(args) || length(args) == 1 || error(
        "Choose one case or 'all'. Run `julia --project=. main.jl --help` for usage.",
    )

    selection = isempty(args) ? "all" : args[1]
    names = selection == "all" ? sort!(collect(keys(CASES))) : [selection]

    for name in names
        haskey(CASES, name) || error("Unknown case '$name'. Use --list to see available cases.")
        run_case(name)
    end

    println()
    println("Completed $(length(names)) case$(length(names) == 1 ? "" : "s").")
    return nothing
end

main()
