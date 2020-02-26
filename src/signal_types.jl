"""
    GNSSSignal

The most general type of
signal type used in `GNSSTools`.
Used for both `GNSSData` and 
`ReplicaSignal` structs.
"""
abstract type GNSSSignal end


"""
    Data

Structure for holding GNSS signal data.
"""
struct GNSSData{T1,T2} <: GNSSSignal
    file_name::String
    f_s::Float64
    f_if::Float64
    t_length::Float64
    start_data_idx::Int64
    t::Array{Float64,1}
    data::Array{Complex{Float64},1}
    data_type::String
    data_start_time::T1
    site_loc_lla::T2
    sample_num::Int64
    total_data_length::Float64
    nADC::Int64
end


"""
    ReplicaSignal

A abstract struct for the replica signal
structs. For use when specifying types
for method arguments.
"""
abstract type ReplicaSignal <: GNSSSignal end



"""
    L1CASignal

Struct for holding L1 C/A GNSS signal properties for
signal generation.
"""
mutable struct L1CASignal <: ReplicaSignal
    type::Val{:l1ca}
    prn::Int64
    f_s::Float64
    t_length::Float64
    f_if::Float64
    f_d::Float64
    fd_rate::Float64
    Tsys::Float64
    CN0::Float64
    ϕ::Float64
    nADC::Int64
    B::Float64
    code_start_idx::Float64
    l1ca_init_code_phase::Float64
    db_init_code_phase::Float64
    t::Array{Float64,1}
    data::Array{Complex{Float64},1}
    include_carrier::Bool
    include_adc::Bool
    include_noise::Bool
    include_databits::Bool
    f_l1ca_d::Float64
    f_l1ca_dd::Float64
    f_db_d::Float64
    f_db_dd::Float64
    sample_num::Int64
    isreplica::Bool
    noexp::Bool
    chipping_rate::Float64
    sig_freq::Float64
    code_length::Int64
    db_code_length::Int64
    databits::Array{Int64,1}
end


"""
    l5q

Array of symbols specifying that a signal
is a an L5 dataless signal.
"""
const l5q = [:l5, :q]


"""
    l5i

Array of symbols specifying that a signal
is a an L5 signal with data bits.
"""
const l5i = [:l5, :i]


"""
    L5QSignal()

Struct for holding L5Q GNSS signal properties for
signal generation.
"""
mutable struct L5QSignal <: ReplicaSignal
    type::Val{:l5q}
    prn::Int64
    f_s::Float64
    t_length::Float64
    f_if::Float64
    f_d::Float64
    fd_rate::Float64
    Tsys::Float64
    CN0::Float64
    ϕ::Float64
    nADC::Int64
    B::Float64
    code_start_idx::Float64
    l5q_init_code_phase::Float64
    nh_init_code_phase::Float64
    t::Array{Float64,1}
    data::Array{Complex{Float64},1}
    include_carrier::Bool
    include_adc::Bool
    include_noise::Bool
    f_l5q_d::Float64
    f_l5q_dd::Float64
    f_nh_d::Float64
    f_nh_dd::Float64
    sample_num::Int64
    isreplica::Bool
    noexp::Bool
    chipping_rate::Float64
    sig_freq::Float64
    code_length::Int64
end


"""
    L5ISignal()

Struct for holding L5Q GNSS signal properties for
signal generation.
"""
mutable struct L5ISignal <: ReplicaSignal
    type::Val{:l5i}
    prn::Int64
    f_s::Float64
    t_length::Float64
    f_if::Float64
    f_d::Float64
    fd_rate::Float64
    Tsys::Float64
    CN0::Float64
    ϕ::Float64
    nADC::Int64
    B::Float64
    code_start_idx::Float64
    l5i_init_code_phase::Float64
    nh_init_code_phase::Float64
    databit_init_code_phase::Float64
    t::Array{Float64,1}
    data::Array{Complex{Float64},1}
    include_carrier::Bool
    include_adc::Bool
    include_noise::Bool
    include_databits::Bool
    f_l5i_d::Float64
    f_l5i_dd::Float64
    f_nh_d::Float64
    f_nh_dd::Float64
    f_db_d::Float64  # databit Doppler
    f_db_dd::Float64
    sample_num::Int64
    isreplica::Bool
    noexp::Bool
    chipping_rate::Float64
    sig_freq::Float64
    code_length::Float64
    db_code_length::Int64
    databits::Array{Int64,1}
end