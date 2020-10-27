"""
    SignalType

Struct for holding parms for a
given signal type.
"""
struct SignalType{T1,T2,T3,T4,T5,T6,T7}
    name::String           # Name of signal type (i.e. "l1ca" or "l5i" or "l5q")
    code_num::T1           # Number of codes in signal type
                           # Includes primary, secondary, and nav codes
    codes::T3              # Vector of dictionaries
                           # Each dictionary is a primary, secondary, or
                           # or nav code which is specific for each PRN number
                           # The last dictionary in the vector of `codes` is
                           # assumed to be the navbits. This "code" will be
                           # ignored when generating replica signals for
                           # acquisition, fine acquisition, and signal tracking
                           # methods. Primary codes are considered the
                           # first code in `codes`.
    chipping_rates::T4     # Vector of chipping rates of length `code_num`
    code_lengths::T5       # Vector of code lengths of length `code_num`
    channel::T6            # Whether signal is on the I or Q channel
                           # If signal is on I, `channel = 1 + 0im`
                           # If signal is on Q, `channel = 0 + 1im`
                           # If signal is on both, `channel = 1 + 1im`
                           # `channel` is multiplied on resulting code value
                           # during signal generation
    include_codes::T7      # Bool vector of length `code_num`
                           # Used to determine if a given code is used for
                           # signal generation.
                           # Will be used in `calc_code_val` method to determine
                           # if a code value will be calculated and XOR`ed to
                           # other code values.
end


"""
    definesignaltype(codes, chipping_rates, code_lengths, channel)

`channel` can be either `"I"`, `"Q"`, or "`both`".
"""
function definesignaltype(codes, chipping_rates, code_lengths, channel="both";
                          name="custom")
    if typeof(codes[1]) ~= Dict
        error("Primary code must be a dictionary of codes. The primary code is the first index of `codes`.")
    end
    code_num = length(codes)
    include_codes = fill(true, code_num)
    if channel == "I"
        channel = 1 + 0im
    elseif channel == "Q"
        channel = 0 + 1im
    elseif channel == "both"
        channel = 1 + 1im
    else
        error("Invalid channel specified.")
    end
    code_keys = collect(keys(codes[1]))
    dict_type = Dict{eltype(code_keys),Vector{eltype(codes[1][code_keys[1]])}}
    signal_codes = Vector{dict_type}(undef, code_num)
    for i in 1:code_num
        if typeof(codes[i]) == Dict
            signal_codes[i] = view(codes, i)
        else
            code = dict_type()
            for code_key in code_keys
                code[code_key] = view(codes, i)
            end
            signal_codes[i] = code
        end
    end
end


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
    ReplicaSignals

A abstract struct for the replica signal
structs. For use when specifying types
for method arguments.
"""
mutable struct ReplicaSignals{T1,T2,T3,T4} <: GNSSSignal
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
    init_code_phases::T1
    t::Array{Float64,1}
    data::Array{Complex{Float64},1}
    include_carrier::Bool
    include_adc::Bool
    include_noise::Bool
    include_databits::Bool
    chipping_rates_d::T2
    chipping_rates_dd::T3
    sample_num::Int64
    isreplica::Bool
    noexp::Bool
    sig_freq::Float64
    signal_type::T4
    # thermal_noise::Array{Complex{Float64},1}
    # phase_noise::Array{Complex{Float64},1}
end



#------------------------------------------------------------------------------
#                             OLD IMPLEMENTATON
#------------------------------------------------------------------------------


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
    # phase_noise::Array{Complex{Float64},1}
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
    # phase_noise::Array{Complex{Float64},1}
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
    db_init_code_phase::Float64
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
    code_length::Int64
    db_code_length::Int64
    databits::Array{Int64,1}
    # phase_noise::Array{Complex{Float64},1}
end
