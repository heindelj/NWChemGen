function parse_nwchem_output_file(output_file::String, observables::String...)
    
end

function get_matching_key(line::String, search_terms::Dict{String, String})
    for key in keys(search_terms)
        if occursin(search_terms[key], line)
            return key
        end
    end
    return ""
end

function get_energies(output_file_lines::Vector{String})
    search_terms = (Dict{String, String}(
        "scf"     => "Total SCF energy",
        "mp2"     => "Total MP2 energy",
        "ccsd(t)" => "Total CCSD(T) energy"
    ))
    energies = Dict{String, Vector{Float64}}()
    for (i, line) in enumerate(output_file_lines)
        possible_key = get_matching_key(line, search_terms)
        if possible_key != ""
            energy = parse(Float64, split(line)[end])

            if haskey(energies, possible_key)
                push!(energies[possible_key], energy)
            else
                energies[possible_key] = [energy]
            end
        end
    end
    return energies
end

function get_gradients(output_file_lines::Vector{String})
    search_terms = (Dict{String, String}(
        "scf"     => "RHF ENERGY GRADIENTS",
        "mp2"     => "mp2 ENERGY GRADIENTS",
        "ccsd(t)" => "CCSD ENERGY GRADIENTS"
    ))
    gradients = Dict{String, Vector{Matrix{Float64}}}()
    for (i, line) in enumerate(output_file_lines)
        possible_key = get_matching_key(line, search_terms)
        if possible_key != ""
            num_atoms::Int = 0
            j::Int = 0
            while output_file_lines[i+j+4] != ""
                num_atoms += 1
                j += 1
            end
            grads = zeros(3, num_atoms)
            j = 0
            while output_file_lines[i+j+4] != ""
                grads[:, j+1] = parse.(Float64, split(nwchem_out[i+j+4])[end-2:end])
                j += 1
            end

            if haskey(gradients, possible_key)
                push!(gradients[possible_key], grads)
            else
                gradients[possible_key] = [grads]
            end
        end
    end
    return gradients
end
