export Fragment, fragments, fragments_df, eachfragment, nfragments, parent_fragment

struct Fragment{T}
    sys::System{T}
    row::DataFrameRow
end

function Fragment(
    chain::Chain{T}, 
    number::Int, 
    name::String = "", 
    properties::Properties = Properties()
) where T
    sys = chain.sys
    idx = _next_idx(sys)
    push!(sys.fragments, (idx, number, name, properties, chain.row.molecule_id, chain.idx))
    _fragment_by_idx(sys, idx)
end

function Base.getproperty(frag::Fragment, name::Symbol)
    in(name, fieldnames(FragmentTuple)) && return getproperty(getfield(frag, :row), name)
    getfield(frag, name)
end

function Base.setproperty!(frag::Fragment, name::Symbol, val)
    in(name, fieldnames(FragmentTuple)) && return setproperty!(getfield(frag, :row), name, val)
    setfield!(frag, name, val)
end

# TODO hide internals
@inline Base.show(io::IO, ::MIME"text/plain", frag::Fragment) = show(io, getfield(frag, :row))
@inline Base.show(io::IO, frag::Fragment) = show(io, getfield(frag, :row))

@inline Base.parent(frag::Fragment) = frag.sys
@inline parent_system(frag::Fragment) = parent(frag)
@inline parent_molecule(frag::Fragment) = _molecule_by_idx(frag.sys, frag.row.molecule_id)
@inline parent_chain(frag::Fragment) = _chain_by_idx(frag.sys, frag.row.chain_id)

@inline function _fragment_by_idx(sys::System{T}, idx::Int) where T
    Fragment{T}(sys, DataFrameRow(sys.fragments, findfirst(sys.fragments.idx .== idx), :))
end

function _fragments(sys::System{T};
        molecule_id::Union{Nothing, Int} = nothing,
        chain_id::Union{Nothing, Int} = nothing
) where T
    isnothing(molecule_id) && isnothing(chain_id) && return sys.fragments

    cols = Tuple{Symbol, Int}[]
    isnothing(molecule_id) || push!(cols, (:molecule_id, molecule_id))
    isnothing(chain_id)    || push!(cols, (:chain_id, chain_id))

    get(
        groupby(sys.fragments, getindex.(cols, 1)),
        ntuple(i -> cols[i][2], length(cols)),
        DataFrame(_SystemFragmentTuple[])
    )
end

@inline function fragments(sys::System; kwargs...)
    collect(eachfragment(sys; kwargs...))
end

@inline function fragments_df(sys::System; kwargs...)
    SystemDataFrame(sys, view(_fragments(sys; kwargs...), :, 1:length(fieldnames(FragmentTuple))))
end

@inline function eachfragment(sys::System{T}; kwargs...) where T
    (Fragment{T}(sys, row) for row in eachrow(_fragments(sys; kwargs...)))
end

@inline function nfragments(sys::System; kwargs...)
    nrow(_fragments(sys; kwargs...))
end

#=
    Molecule fragments
=#
@inline _fragments(mol::Molecule; kwargs...) = _fragments(mol.sys; molecule_id = mol.idx, kwargs...)
@inline fragments(mol::Molecule; kwargs...) = fragments(mol.sys; molecule_id = mol.idx, kwargs...)
@inline fragments_df(mol::Molecule; kwargs...) = fragments_df(mol.sys; molecule_id = mol.idx, kwargs...)
@inline eachfragment(mol::Molecule; kwargs...) = eachfragment(mol.sys; molecule_id = mol.idx, kwargs...)
@inline nfragments(mol::Molecule; kwargs...) = nfragments(mol.sys; molecule_id = mol.idx, kwargs...)

#=
    Chain fragments
=#
@inline _fragments(chain::Chain; kwargs...) = _fragments(chain.sys; chain_id = chain.idx, kwargs...)
@inline fragments(chain::Chain; kwargs...) = fragments(chain.sys; chain_id = chain.idx, kwargs...)
@inline fragments_df(chain::Chain; kwargs...) = fragments_df(chain.sys; chain_id = chain.idx, kwargs...)
@inline eachfragment(chain::Chain; kwargs...) = eachfragment(chain.sys; chain_id = chain.idx, kwargs...)
@inline nfragments(chain::Chain; kwargs...) = nfragments(chain.sys; chain_id = chain.idx, kwargs...)

@inline function Base.push!(chain::Chain, frag::FragmentTuple)
    push!(chain.sys.fragments, (_with_idx(frag, _next_idx(chain.sys))..., chain.idx))
    chain
end

#=
    Fragment atoms
=#
@inline _atoms(frag::Fragment; kwargs...) = _atoms(frag.sys; fragment_id = frag.idx, kwargs...)
@inline atoms(frag::Fragment; kwargs...) = atoms(frag.sys; fragment_id = frag.idx, kwargs...)
@inline atoms_df(frag::Fragment; kwargs...) = atoms_df(frag.sys; fragment_id = frag.idx, kwargs...)
@inline eachatom(frag::Fragment; kwargs...) = eachatom(frag.sys; fragment_id = frag.idx, kwargs...)
@inline natoms(frag::Fragment; kwargs...) = natoms(frag.sys; fragment_id = frag.idx, kwargs...)

@inline function Base.push!(frag::Fragment{T}, atom::AtomTuple{T}; kwargs...) where T
    push!(frag.sys, atom; molecule_id = frag.row.molecule_id, chain_id = frag.row.chain_id,
        fragment_id = frag.idx, kwargs...)
    frag
end

#=
    Fragment bonds
=#
@inline _bonds(frag::Fragment; kwargs...) = _bonds(frag.sys; fragment_id = frag.idx, kwargs...)
@inline bonds(frag::Fragment; kwargs...) = bonds(frag.sys; fragment_id = frag.idx, kwargs...)
@inline bonds_df(frag::Fragment; kwargs...) = bonds_df(frag.sys; fragment_id = frag.idx, kwargs...)
@inline eachbond(frag::Fragment; kwargs...) = eachbond(frag.sys; fragment_id = frag.idx, kwargs...)
@inline nbonds(frag::Fragment; kwargs...) = nbonds(frag.sys; fragment_id = frag.idx, kwargs...)

@inline function Base.push!(frag::Fragment, bond::Bond)
    push!(frag.sys, bond)
    frag
end
