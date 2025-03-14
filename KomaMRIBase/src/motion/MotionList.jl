include("Interpolation.jl")
include("SpinSpan.jl")
include("TimeCurve.jl")
include("Action.jl")
include("Motion.jl")

"""
    motionlist = MotionList(motions...)

MotionList struct. The other option, instead of `NoMotion`, 
is to define a dynamic phantom by means of the `MotionList` struct.
It is composed by one or more [`Motion`](@ref) instances. 

# Arguments
- `motions`: (`::Vector{Motion{T<:Real}}`) vector of `Motion` instances

# Returns
- `motionlist`: (`::MotionList`) MotionList struct

# Examples
```julia-repl
julia>  motionlist = MotionList(
            Motion(
                action = Translate(0.01, 0.0, 0.02),
                time = TimeRange(0.0, 1.0),
                spins = AllSpins()
            ),
            Motion(
                action = Rotate(0.0, 0.0, 45.0),
                time = Periodic(1.0),
                spins = SpinRange(1:10)
            )
        )
```
"""
struct MotionList{T<:Real}
    motions::Vector{<:Motion{T}}
end

# NOTE: this constructor must be simplified once the Vector{<:Motion} approach is accomplished: 
# https://github.com/JuliaHealth/KomaMRI.jl/issues/480
""" Constructors """
function MotionList(motions::Motion...) 
    if length(motions) == 0
        return NoMotion()
    elseif length(motions) == 1
        return motions[1]
    else
        return MotionList([motions...])
    end 
end

# NOTE: these vcat methods must be simplified once the Vector{<:Motion} approach is accomplished: 
# https://github.com/JuliaHealth/KomaMRI.jl/issues/480
""" Addition of MotionLists """ 
# MotionList + MotionList
function Base.vcat(m1::MotionList{T}, m2::MotionList{T}, Ns1, Ns2) where {T<:Real}
    mv_aux = Motion{T}[]
    for m in m1.motions
        m_aux = copy(m)
        m_aux.spins = expand(m_aux.spins, Ns1)
        push!(mv_aux, m_aux)
    end
    for m in m2.motions
        m_aux = copy(m)
        m_aux.spins = expand(m_aux.spins, Ns2)
        m_aux.spins = SpinRange(m_aux.spins.range .+ Ns1)
        push!(mv_aux, m_aux)
    end
    return MotionList(mv_aux...)
end
# Motion + Motion
function Base.vcat(m1::Motion{T}, m2::Motion{T}, Ns1, Ns2) where {T<:Real}
    mv_aux = Motion{T}[]
    m_aux = copy(m1)
    m_aux.spins = expand(m_aux.spins, Ns1)
    push!(mv_aux, m_aux)
    m_aux = copy(m2)
    m_aux.spins = expand(m_aux.spins, Ns2)
    m_aux.spins = SpinRange(m_aux.spins.range .+ Ns1)
    push!(mv_aux, m_aux)
    return MotionList(mv_aux...)
end
# Motion + MotionList
Base.vcat(m1::MotionList{T}, m2::Motion{T}, Ns1, Ns2) where {T<:Real} = vcat(m2, m1, Ns2, Ns1)
function Base.vcat(m1::Motion{T}, m2::MotionList{T}, Ns1, Ns2) where {T<:Real}
    mv_aux = Motion{T}[]
    m_aux = copy(m1)
    m_aux.spins = expand(m_aux.spins, Ns1)
    push!(mv_aux, m_aux)
    for m in m2.motions
        m_aux = copy(m)
        m_aux.spins = expand(m_aux.spins, Ns2)
        m_aux.spins = SpinRange(m_aux.spins.range .+ Ns1)
        push!(mv_aux, m_aux)
    end
    return MotionList(mv_aux...)
end

""" MotionList sub-group """
function Base.getindex(mv::MotionList{T}, p) where {T<:Real}
    motion_array_aux = Motion{T}[]
    for m in mv.motions
        m[p] isa NoMotion ? nothing : push!(motion_array_aux, m[p])
    end
    return MotionList(motion_array_aux...)
end
function Base.view(mv::MotionList{T}, p) where {T<:Real}
    motion_array_aux = Motion{T}[]
    for m in mv.motions
        @view(m[p]) isa NoMotion ? nothing : push!(motion_array_aux, @view(m[p]))
    end
    return MotionList(motion_array_aux...)
end

""" Compare two MotionLists """
function Base.:(==)(mv1::MotionList{T}, mv2::MotionList{T}) where {T<:Real}
    if length(mv1) != length(mv2) return false end
    sort_motions!(mv1)
    sort_motions!(mv2)
    return reduce(&, mv1.motions .== mv2.motions)
end
function Base.:(≈)(mv1::MotionList{T}, mv2::MotionList{T}) where {T<:Real} 
    if length(mv1) != length(mv2) return false end
    sort_motions!(mv1)
    sort_motions!(mv2)
    return reduce(&, mv1.motions .≈ mv2.motions)
end

""" MotionList length """
Base.length(m::MotionList) = length(m.motions)

"""
    x, y, z = get_spin_coords(motionset, x, y, z, t)

Calculates the position of each spin at a set of arbitrary time instants, i.e. the time steps of the simulation. 
For each dimension (x, y, z), the output matrix has ``N_{\t{spins}}`` rows and `length(t)` columns.

# Arguments
- `motion`: (`::Union{NoMotion, MotionList{T<:Real}}`) phantom motion
- `x`: (`::AbstractVector{T<:Real}`, `[m]`) spin x-position vector
- `y`: (`::AbstractVector{T<:Real}`, `[m]`) spin y-position vector
- `z`: (`::AbstractVector{T<:Real}`, `[m]`) spin z-position vector
- `t`: horizontal array of time instants

# Returns
- `x, y, z`: (`::Tuple{AbstractArray, AbstractArray, AbstractArray}`) spin positions over time
"""
function get_spin_coords(
    ml::MotionList{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t
) where {T<:Real}
    # Sort motions
    sort_motions!(ml)
    # Buffers for positions:
    xt, yt, zt = x .+ 0*t, y .+ 0*t, z .+ 0*t
    # Buffers for displacements:
    ux, uy, uz = xt .* zero(T), yt .* zero(T), zt .* zero(T)
    # Composable motions: they need to be run sequentially. Note that they depend on xt, yt, and zt
    for m in Iterators.filter(is_composable, ml.motions)
        t_unit = unit_time(t, m.time.t, m.time.t_unit, m.time.periodic, m.time.periods)
        idx = get_indexing_range(m.spins)
        displacement_x!(@view(ux[idx, :]), m.action, @view(xt[idx, :]), @view(yt[idx, :]), @view(zt[idx, :]), t_unit)
        displacement_y!(@view(uy[idx, :]), m.action, @view(xt[idx, :]), @view(yt[idx, :]), @view(zt[idx, :]), t_unit)
        displacement_z!(@view(uz[idx, :]), m.action, @view(xt[idx, :]), @view(yt[idx, :]), @view(zt[idx, :]), t_unit)
        xt .+= ux; yt .+= uy; zt .+= uz
        ux .*= zero(T); uy .*= zero(T); uz .*= zero(T)
    end
    # Additive motions: these motions can be run in parallel
    for m in Iterators.filter(!is_composable, ml.motions)
        t_unit = unit_time(t, m.time.t, m.time.t_unit, m.time.periodic, m.time.periods)
        idx = get_indexing_range(m.spins)
        displacement_x!(@view(ux[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
        displacement_y!(@view(uy[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
        displacement_z!(@view(uz[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
        xt .+= ux; yt .+= uy; zt .+= uz
        ux .*= zero(T); uy .*= zero(T); uz .*= zero(T)
    end
    return xt, yt, zt
end

"""
    times = times(motion)
"""
function times(ml::MotionList)
    nodes = reduce(vcat, [times(m) for m in ml.motions])
    return unique(sort(nodes))
end

"""
    sort_motions!(motion)

Sorts motions in a list according to their starting time. It modifies the original list.
If `motionset::NoMotion`, this function does nothing.
If `motionset::MotionList`, this function sorts its motions.

# Arguments
- `motion`: (`::Union{NoMotion, MotionList{T<:Real}}`) phantom motion

# Returns
- `nothing`
"""
function sort_motions!(m::MotionList)
    sort!(m.motions; by=m -> m.time.t_start)
    return nothing
end

function add_jump_times!(t, ml::MotionList)
    for m in ml.motions
        add_jump_times!(t, m)
    end
end