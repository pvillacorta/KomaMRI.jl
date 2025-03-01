"""
    motion = Motion(action)
    motion = Motion(action, time)
    motion = Motion(action, time, spins)

Motion struct. It defines the motion, during a certain time interval,
of a given group of spins. It is composed by three fields: `action`, which 
defines the motion itself, `time`, which accounts for the time during
which the motion takes place, and `spins`, which indicates the spins 
that are affected by that motion.

# Arguments
- `action`: (`::AbstractAction{T<:Real}`) action, such as [`Translate`](@ref) or [`Rotate`](@ref)
- `time`: (`::TimeCurve{T<:Real}`, `=TimeRange(0.0)`) time information about the motion
- `spins`: (`::AbstractSpinSpan`, `=AllSpins()`) spin indexes affected by the motion

# Returns
- `motion`: (`::Motion`) Motion struct

# Examples
```julia-repl
julia> motion =  Motion(
            action = Translate(0.01, 0.0, 0.02),
            time = TimeRange(0.0, 1.0),
            spins = SpinRange(1:10)
       )
```
"""
@with_kw mutable struct Motion{T<:Real}
    action::AbstractAction{T}
    time  ::TimeCurve{T}      = TimeRange(t_start=zero(typeof(action).parameters[1]), t_end=eps(typeof(action).parameters[1]))
    spins ::AbstractSpinSpan  = AllSpins()
end

# Main constructors
function Motion(action) 
    T = first(typeof(action).parameters)
    return Motion(action, TimeRange(t_start=zero(T), t_end=eps(T)), AllSpins())
end
function Motion(action, time::TimeCurve)
    T = first(typeof(action).parameters)
    return Motion(action, time, AllSpins())
end
function Motion(action, spins::AbstractSpinSpan)
    T = first(typeof(action).parameters)
    return Motion(action, TimeRange(t_start=zero(T), t_end=eps(T)), spins)
end

# Custom constructors
"""
    translate = Translate(dx, dy, dz, time, spins)

# Arguments
- `dx`: (`::Real`, `[m]`) translation in x
- `dy`: (`::Real`, `[m]`) translation in y 
- `dz`: (`::Real`, `[m]`) translation in z
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `translate`: (`::Motion`) Motion struct

# Examples
```julia-repl
julia> translate = Translate(0.01, 0.02, 0.03, TimeRange(0.0, 1.0), SpinRange(1:10))
```
"""
function Translate(dx, dy, dz, time=TimeRange(t_start=zero(eltype(dx)), t_end=eps(eltype(dx))), spins=AllSpins())
    return Motion(Translate(dx, dy, dz), time, spins)
end

"""
    rotate = Rotate(pitch, roll, yaw, spins)

# Arguments
- `pitch`: (`::Real`, `[º]`) rotation in x
- `roll`: (`::Real`, `[º]`) rotation in y 
- `yaw`: (`::Real`, `[º]`) rotation in z
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `rotate`: (`::Motion`) Motion struct with [`Rotate`](@ref) action

# Examples
```julia-repl
julia> rotate = Rotate(15.0, 0.0, 20.0, TimeRange(0.0, 1.0), SpinRange(1:10))
```
"""
function Rotate(pitch, roll, yaw, time=TimeRange(t_start=zero(eltype(pitch)), t_end=eps(eltype(pitch))), spins=AllSpins())
    return Motion(Rotate(pitch, roll, yaw), time, spins)
end

"""
    heartbeat = HeartBeat(circumferential_strain, radial_strain, longitudinal_strainl, time, spins)

# Arguments
- `circumferential_strain`: (`::Real`) contraction parameter
- `radial_strain`: (`::Real`) contraction parameter
- `longitudinal_strain`: (`::Real`) contraction parameter
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `heartbeat`: (`::Motion`) Motion struct with [`HeartBeat`](@ref) action

# Examples
```julia-repl
julia> heartbeat = HeartBeat(-0.3, -0.2, 0.0, TimeRange(0.0, 1.0), SpinRange(1:10))
```
"""
function HeartBeat(circumferential_strain, radial_strain, longitudinal_strain, time=TimeRange(t_start=zero(eltype(circumferential_strain)), t_end=eps(eltype(circumferential_strain))), spins=AllSpins())
    return Motion(HeartBeat(circumferential_strain, radial_strain, longitudinal_strain), time, spins)
end

"""
    path = Path(dx, dy, dz, time, spins)

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `path`: (`::Motion`) Motion struct with [`Path`](@ref) action

# Examples
```julia-repl
julia> path = Path(
          [0.01 0.02; 0.02 0.03], 
          [0.02 0.03; 0.03 0.04], 
          [0.03 0.04; 0.04 0.05], 
          TimeRange(0.0, 1.0), 
          SpinRange(1:10)
       )
```
"""
function Path(dx, dy, dz, time=TimeRange(t_start=zero(eltype(dx)), t_end=eps(eltype(dx))), spins=AllSpins())
    return Motion(Path(dx, dy, dz), time, spins)
end

"""
    flowpath = FlowPath(dx, dy, dz, spin_reset, time, spins)

# Arguments
- `dx`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in x
- `dy`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in y 
- `dz`: (`::AbstractArray{T<:Real}`, `[m]`) displacements in z
- `spin_reset`: (`::AbstractArray{Bool}`) reset spin state flags
- `time`: (`::TimeCurve{T<:Real}`) time information about the motion
- `spins`: (`::AbstractSpinSpan`) spin indexes affected by the motion

# Returns
- `flowpath`: (`::Motion`) Motion struct with [`FlowPath`](@ref) action

# Examples
```julia-repl
julia> flowpath = FlowPath(
          [0.01 0.02; 0.02 0.03], 
          [0.02 0.03; 0.03 0.04], 
          [0.03 0.04; 0.04 0.05], 
          [false false; false true],
          TimeRange(0.0, 1.0), 
          SpinRange(1:10)
       )
```
"""
function FlowPath(dx, dy, dz, spin_reset, time=TimeRange(t_start=zero(eltype(dx)), t_end=eps(eltype(dx))), spins=AllSpins())
    return Motion(FlowPath(dx, dy, dz, spin_reset), time, spins)
end

""" Compare two Motions """
Base.:(==)(m1::Motion, m2::Motion) = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field) == getfield(m2, field) for field in fieldnames(typeof(m1))])
Base.:(≈)(m1::Motion, m2::Motion)  = (typeof(m1) == typeof(m2)) & reduce(&, [getfield(m1, field)  ≈ getfield(m2, field) for field in fieldnames(typeof(m1))])

""" Motion sub-group """
function Base.getindex(m::Motion, p)
    idx, spin_range = m.spins[p]
    return idx !== nothing ? Motion(m.action[idx], m.time, spin_range) : NoMotion()
end
function Base.view(m::Motion, p)
    idx, spin_range = @view(m.spins[p])
    return idx !== nothing ? Motion(@view(m.action[idx]), m.time, spin_range) : NoMotion()
end

"""
    x, y, z = get_spin_coords(motion, x, y, z, t)
"""
function get_spin_coords(
    m::Motion{T}, x::AbstractVector{T}, y::AbstractVector{T}, z::AbstractVector{T}, t
) where {T<:Real}
    ux, uy, uz = x .* (0*t), y .* (0*t), z .* (0*t) # Buffers for displacements
    t_unit = unit_time(t, m.time.t, m.time.t_unit, m.time.periodic, m.time.periods)
    idx = get_indexing_range(m.spins)
    displacement_x!(@view(ux[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
    displacement_y!(@view(uy[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
    displacement_z!(@view(uz[idx, :]), m.action, @view(x[idx]), @view(y[idx]), @view(z[idx]), t_unit)
    return x .+ ux, y .+ uy, z .+ uz
end

# Auxiliary functions
times(m::Motion) = times(m.time.t, m.time.periods)
is_composable(m::Motion) = is_composable(m.action)
add_jump_times!(t, m::Motion) = add_jump_times!(t, m.action, m.time)
add_jump_times!(t, ::AbstractAction, ::TimeCurve) = nothing