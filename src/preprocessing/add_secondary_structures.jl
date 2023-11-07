export add_secondary_structures!


function add_secondary_structures!(system::System{T}, pdb_path::String) where {T<:Real}
    structure_list_per_chain = Dict()
    for chain in eachchain(system)
        structure_list_per_chain[chain.name] = []
    end
    found_section = false
    file = open(pdb_path)
    for line in eachline(file)
        if(startswith(line, "HELIX"))
            found_section = true
            chainname = line[20:20]
            start_number = tryparse(Int, line[22:25])
            end_number = tryparse(Int, line[34:37])
            if(start_number!==nothing && end_number!==nothing)
                push!(structure_list_per_chain[chainname], (start_number, end_number, SecondaryStructure.HELIX))
            end
        elseif(startswith(line, "SHEET"))
            found_section = true
            chainname = line[22:22]
            start_number = tryparse(Int, line[23:26])
            end_number = tryparse(Int, line[34:37])
            if(start_number!==nothing && end_number!==nothing)
                push!(structure_list_per_chain[chainname], (start_number, end_number, SecondaryStructure.SHEET))
            end
        else # stop reading the file when the secondary structure section of the pdb file is finished
            if(found_section) 
                break
            end
        end
    end
    close(file)


    for chain in eachchain(system)
        for fragment in eachfragment(chain)
            assigned_structure = false
            for (start_index, end_index, structure) in structure_list_per_chain[chain.name]
                if(start_index<=fragment.number && fragment.number<=end_index)
                    fragment.properties[:SS] = structure
                    assigned_structure = true
                    break
                end
            end
            if(!assigned_structure)
                fragment.properties[:SS] = SecondaryStructure.NONE
            end
        end
        chain.properties[:SSs] = structure_list_per_chain[chain.name]
    end
end

function add_secondary_structures!(system::System{T}, orig_pdb::ProteinStructure) where {T<:Real} 
    #TODO write tests, comments, think about different types of helices
    structure_mapping = Dict(
        'H'=>SecondaryStructure.HELIX, 
        'B'=>SecondaryStructure.SHEET, 
        'E'=>SecondaryStructure.SHEET, 
        'G'=>SecondaryStructure.HELIX, 
        'I'=>SecondaryStructure.HELIX, 
        'P'=>SecondaryStructure.HELIX, 
        'T'=>SecondaryStructure.NONE, 
        'S'=>SecondaryStructure.NONE, 
        ' '=>SecondaryStructure.NONE, 
        '-'=>SecondaryStructure.NONE)

    for chain in eachchain(system)
        bs_chain = orig_pdb[chain.name]
        structure_list = []
        current_structure = nothing
        start_index = -1

        for residue in bs_chain
            if(sscode(residue)!=current_structure)
                # end last ss
                if(current_structure!==nothing 
                    && start_index>0 
                    && structure_mapping[current_structure]!=SecondaryStructure.NONE)
                    push!(structure_list, (start_index, resnumber(residue)-1, structure_mapping[current_structure]))
                end

                # start new ss
                current_structure = sscode(residue)
                start_index = resnumber(residue)
            end
        end
        # end last structure
        if(current_structure!==nothing 
            && start_index>0 
            && structure_mapping[current_structure]!=SecondaryStructure.NONE)
            last_residue = chain[lastindex(chain)]
            push!(structure_list, (start_index, resnumber(last_residue), structure_mapping[current_structure]))
        end

        for fragment in eachfragment(chain)
            assigned_structure = false
            for (start_index, end_index, structure) in structure_list
                if(start_index<=fragment.number && fragment.number<=end_index)
                    fragment.properties[:SS] = structure
                    assigned_structure = true
                end
            end
            if(!assigned_structure)
                fragment.properties[:SS] = SecondaryStructure.NONE
            end
        end
        chain.properties[:SSs] = structure_list
    end
    return system
end