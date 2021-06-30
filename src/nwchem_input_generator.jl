
Base.@kwdef mutable struct NWChemInput
    basis::Dict{String, String} # for specifying different basis sets for different atoms
    theory::String
    task::Vector{String} = ["energy"] # can specify multiple tasks
    memory::Int = 1000 # mb
    block_settings::Vector{Tuple{String, Dict{String, String}}} = []
    special_settings::Vector{String} = [] # for e.g. set n_lin_dep 0
    geometry_block_settings::Vector{String} = ["noautosym", "noautoz"]
    header_string::String    = ""           # the string for the header info like memory, etc.
    geometry_strings::Vector{String} = []   # the string for the geometry block
    settings_string::String  = ""           # the string for everything after the geometry.
end

NWChemInput(basis::Dict{String, String}, theory::String) = NWChemInput(basis = basis, theory = theory)
NWChemInput(basis::String, theory::String) = NWChemInput(basis = Dict("*" => basis), theory = theory)

function write_input_file(input::NWChemInput, geoms::Vector{Matrix{Float64}}, atom_labels::Vector{Vector{String}}, input_file_name::String="input.nw", out_directory::String="nwchem")
    set_header_string!(input)
    set_geometry_string!(input, geoms, atom_labels)
    set_settings_string!(input)
    input_file::String = input.header_string 
    for geom_string in input.geometry_strings
        input_file = string(input_file, geom_string, input.settings_string)
    end
    
    if !isdir(out_directory)
        mkdir(out_directory)
    end
    open(next_unique_name(string(out_directory, "/", input_file_name)), "w") do io
        write(io, input_file)
    end
end

function next_unique_name(file_name::String, i::Int=1)
    if isfile(string(file_name))
        new_file_name = string(splitext(file_name)[1], "_", i, splitext(file_name)[2])
        if isfile(new_file_name)
            new_file_name = next_unique_name(file_name, i + 1)
        end
        return new_file_name
    end
    return file_name
end

function geometry_to_string(geom::AbstractMatrix{Float64}, atom_labels::Vector{String})
    geom_string = ""
    for (i, vec) in enumerate(eachcol(geom))
        geom_string = string(geom_string, atom_labels[i], " ", join(vec, " "), "\n")
    end
    return geom_string
end

function set_geometry_string!(input::NWChemInput, geoms::Vector{Matrix{Float64}}, atom_labels::Vector{Vector{String}})
    input.geometry_strings = [string("geometry ", join(input.geometry_block_settings, " "), "\n", geometry_to_string(geoms[i], atom_labels[i]), "end\n") for i in 1:length(geoms)]
end

function set_settings_string!(input::NWChemInput)
    input.settings_string = string(get_basis_string(input), get_special_settings(input), get_block_strings(input), get_task_string(input))
end

function set_header_string!(input::NWChemInput)
    input.header_string = string("echo\n", get_memory_string(input))
end

function set_basis!(input::NWChemInput, basis::Union{String, Dict{String, String}})
    if typeof(basis) == String
        input.basis = Dict("*" => basis)
    else
        input.basis = Dict([(lowercase(key), lowercase(val)) for (key, val) in
pairs(basis)])
    end
end

function get_basis_string(input::NWChemInput)
    return string("basis spherical\n ", join([string(atom, " library ", basis, "\n") for (atom, basis) in pairs(input.basis)]), "end\n")
end

function set_memory!(input::NWChemInput, memory::Int)
    @assert memory > 0 "Need a non-zero amount of memory in MB."
    input.memory = memory
end

function get_memory_string(input::NWChemInput)
    return string("memory ", input.memory, " mb\n")
end

function set_theory!(input::NWChemInput, theory::String)
    possible_theories = ["hf", "scf", "mp2", "ccsd", "ccsd(t)", "dft"]
    theory = lowercase(theory)
    if theory in possible_theories
        theory == "hf" ? input.theory = "scf" : input.theory = theory
    else
        println(string("WARNING: Requested theoretical method ", theory, " is not known. Proceeding, but check that this is valid input before running."))
        input.theory = theory
    end
end

function set_task!(input::NWChemInput, task::Union{String, Vector{String}})
    append!(input.task, task)
end

function get_task_string(input::NWChemInput)
    return join([string("task ", input.theory, " ", task, "\n") for task in input.task])
end

function set_block!(input::NWChemInput, block_name::String, block_settings::Pair{String, String}...)
    settings_dict = Dict{String, String}()
    for setting in block_settings
        push!(settings_dict, setting)
    end
    input.block_settings = (block_name, settings_dict)
end

function get_block_string(block::Tuple{String, Dict{String, String}})
    return string(block[0], "\n", join([string(key, " ", val, "\n") for (key, val) in pairs(input.basis)]), "end\n")
end

function get_block_strings(input::NWChemInput)
    if isempty(input.block_settings)
        return ""
    end
    return join(get_block_string.(input.block_settings))
end

function set_special_settings!(input::NWChemInput, setting::Union{String, Vector{String}})
    append!(input.special_settings, setting)
end

function get_special_settings(input::NWChemInput)
    if isempty(input.special_settings)
        return ""
    end
    return string(join(input.special_settings, "\n"), "\n")
end
