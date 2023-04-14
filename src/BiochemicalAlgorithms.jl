module BiochemicalAlgorithms

using AutoHashEquals
using DataFrames
using DataFramesMeta
using DocStringExtensions

include("core/types.jl")
include("core/element.jl")
include("core/amino_acid.jl")
include("core/bond_order.jl")
include("core/tuples.jl")
include("core/system_component.jl")
include("core/atom_container.jl")
include("core/system.jl")
include("core/atom.jl")
include("core/bond.jl")
include("core/molecule.jl")
include("core/chain.jl")
include("core/fragment.jl")
include("core/nucleotide.jl")
include("core/residue.jl")
include("core/protein.jl")
include("core/moleculargraph_wrapper.jl")

include("substructures/substructure.jl")
include("substructures/smarts.jl")
include("substructures/sssr.jl")

module PubChem
include("fileformats/pubchem_json.jl")
end
include("fileformats/pdb.jl")

include("mappings/atom_bijection.jl")
include("mappings/rigid_mapping.jl")

include("forcefields/MMFF94/mmff94_parameters.jl")

include("preprocessing/fragmentDB.jl")
include("preprocessing/normalize_names.jl")
include("preprocessing/build_bonds.jl")

using .PubChem

export load_pubchem_json, ball_data_path

ball_data_path(parts...) = normpath(joinpath(@__DIR__, "..", "data", parts...))

end
