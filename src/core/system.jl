using AutoHashEquals
export System, SystemDataFrame, default_system, parent_system

"""
    const _SystemAtomTuple{T} = NamedTuple{...}

System-specific extension of `AtomTuple{T}`. See [`AtomTuple`](@ref).

# Additional fields
 - `frame_id::Int`
 - `molecule_id::MaybeInt`
 - `chain_id::MaybeInt`
 - `fragment_id::MaybeInt`
 - `nucleotide_id::MaybeInt`
 - `residue_id::MaybeInt`
"""
const _SystemAtomTuple{T} = NamedTuple{
    (fieldnames(AtomTuple{T})...,
        :frame_id, :molecule_id, :chain_id, :fragment_id, :nucleotide_id, :residue_id),
    Tuple{fieldtypes(AtomTuple{T})...,
        Int, MaybeInt, MaybeInt, MaybeInt, MaybeInt, MaybeInt}
}

"""
    const _SystemChainTuple{T} = NamedTuple{...}

System-specific extension of `ChainTuple{T}`. See [`ChainTuple`](@ref).

# Additional fields
 - `molecule_id::Int`
"""
const _SystemChainTuple = NamedTuple{
    (fieldnames(ChainTuple)..., :molecule_id),
    Tuple{fieldtypes(ChainTuple)..., Int}
}

"""
    const _SystemFragmentTuple{T} = NamedTuple{...}

System-specific extension of `FragmentTuple{T}`. See [`FragmentTuple`](@ref).

# Additional fields
 - `molecule_id::Int`
 - `chain_id::Int`
"""
const _SystemFragmentTuple = NamedTuple{
    (fieldnames(FragmentTuple)..., :molecule_id, :chain_id),
    Tuple{fieldtypes(FragmentTuple)..., Int, Int}
}

"""
    const _SystemNucleotideTuple{T} = NamedTuple{...}

System-specific extension of `NucleotideTuple{T}`. See [`NucleotideTuple`](@ref).

# Additional fields
 - `molecule_id::Int`
 - `chain_id::Int`
"""
const _SystemNucleotideTuple = NamedTuple{
    (fieldnames(NucleotideTuple)..., :molecule_id, :chain_id),
    Tuple{fieldtypes(NucleotideTuple)..., Int, Int}
}

"""
    const _SystemResidueTuple{T} = NamedTuple{...}

System-specific extension of `ResidueTuple{T}`. See [`ResidueTuple`](@ref).

# Additional fields
 - `molecule_id::Int`
 - `chain_id::Int`
"""
const _SystemResidueTuple = NamedTuple{
    (fieldnames(ResidueTuple)..., :molecule_id, :chain_id),
    Tuple{fieldtypes(ResidueTuple)..., Int, Int}
}

"""
    $(TYPEDEF)

Mutable representation of a biomolecular system.

# Fields
 - `name::String`
 - `properties::Properties`

# Constructors
    System(name::String = "", properties::Properties = Properties())

Creates a new and empty `System{Float32}`.

    System{T}(name::String = "", properties::Properties = Properties())

Creates a new and empty `System{T}`.
"""
@auto_hash_equals mutable struct System{T} <: AbstractAtomContainer{T}
    name::String
    properties::Properties

    _atoms::DataFrame
    _bonds::DataFrame
    _molecules::DataFrame
    _chains::DataFrame
    _fragments::DataFrame
    _nucleotides::DataFrame
    _residues::DataFrame
    _curr_idx::Int

    function System{T}(name::String = "", properties::Properties = Properties()) where T
        new(
            name,
            properties,
            DataFrame(_SystemAtomTuple{T}[]),
            DataFrame(BondTuple[]),
            DataFrame(MoleculeTuple[]),
            DataFrame(_SystemChainTuple[]),
            DataFrame(_SystemFragmentTuple[]),
            DataFrame(_SystemNucleotideTuple[]),
            DataFrame(_SystemResidueTuple[]),
            0
        )
    end
end

System(name::String = "", properties::Properties = Properties()) = System{Float32}(name, properties)

"""
    const _default_system

Global default system.
"""
const _default_system = System("default")

"""
    $(TYPEDSIGNATURES)

Returns the global default system.
"""
@inline function default_system()
    _default_system
end

"""
    $(TYPEDSIGNATURES)

Returns the next available `idx` for the given system.
"""
@inline function _next_idx(sys::System{T}) where T
    sys._curr_idx += 1
end

Base.show(io::IO, ::MIME"text/plain", sys::System) = show(io, sys)
Base.show(io::IO, sys::System) = print(io, 
    "System with ", natoms(sys), " atoms", isempty(sys.name) ? "" : " ($(sys.name))")

@doc raw"""
    parent(::Atom)
    parent(::Bond)
    parent(::Chain)
    parent(::Fragment)
    parent(::Molecule)
    parent(::Nucleotide)
    parent(::Protein)
    parent(::Residue)

Returns the `System{T}` containing the given object.
""" Base.parent(::System)

@doc raw"""
    parent_system(::Atom)
    parent_system(::Bond)
    parent_system(::Chain)
    parent_system(::Fragment)
    parent_system(::Molecule)
    parent_system(::Nucleotide)
    parent_system(::Protein)
    parent_system(::Residue)

Returns the `System{T}` containing the given object. Alias for 
[`Base.parent`](@ref Base.parent(::System)).
""" parent_system

"""
    $(TYPEDEF)

System-aware `DataFrame` wrapper.
"""
struct SystemDataFrame{T} <: AbstractDataFrame
    sytem::System{T}
    df::SubDataFrame
end

DataFrames.describe(sys::SystemDataFrame) = describe(getfield(sys, :df))
DataFrames.index(sys::SystemDataFrame) = DataFrames.index(getfield(sys, :df))
DataFrames.nrow(sys::SystemDataFrame) = nrow(getfield(sys, :df))
DataFrames._check_consistency(sys::SystemDataFrame) = DataFrames._check_consistency(getfield(sys, :df))
Base.getindex(sys::SystemDataFrame, row, col) = getindex(getfield(sys, :df), row, col)
Base.getindex(sys::SystemDataFrame, row::Integer, col::Colon) = getindex(getfield(sys, :df), row, col)
Base.setindex!(sys::SystemDataFrame, val, idx) = setindex!(getfield(sys, :df), val, idx)
Base.setindex!(sys::SystemDataFrame, val, row, col) = setindex!(getfield(sys, :df), val, row, col)
Base.getproperty(sys::SystemDataFrame, name::Symbol) = getproperty(getfield(sys, :df), name)
Base.getproperty(sys::SystemDataFrame, name::Symbol, order::Symbol) = getproperty(getfield(sys, :df), name, order)
Base.setproperty!(sys::SystemDataFrame, col, v) = setproperty!(getfield(sys, :df), col, v)
Base.show(io::IO, ::MIME"text/plain", sys::SystemDataFrame) = show(io, getfield(sys, :df))
Base.show(io::IO, sys::SystemDataFrame) = show(io, getfield(sys, :df))
