using BiochemicalAlgorithms
using Test
using DataFrames

@testset verbose=true "BiochemicalAlgorithms.jl" begin
   
    @testset "Core" begin 
        include("test_amino_acid.jl")
        include("test_bond_order.jl")
        include("test_atom.jl")
        include("test_bond.jl")
        include("test_element.jl")
        include("test_fragment.jl")
        include("test_molecule.jl")
        include("test_nucleotide.jl")
        include("test_protein.jl")
        include("test_residue.jl")
        include("test_types.jl")
    end

    @testset "Fileformats" begin 
        include("fileformats/test_pdb.jl")
        include("fileformats/test_pubchem_json.jl")
        include("fileformats/test_sdfile.jl")
    end
    
    @testset "Preprocessing" begin
        include("preprocessing/test_fragmentdb.jl")
    end

    @testset "Substructures" begin
        include("test_substructure.jl")
    end
end
