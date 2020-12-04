"""
    CodeType

Struct for holding parms for a custum code type for either the `I` or `Q`
channel

Fields:

- `name::String`: name of signal type (default = channel (`I`, `Q`, or `IQ`))
- `code_num::Int`: number of codes, including databits, on channel
- `codes::Array{Dict,1}`: array of dictionaries, where each dictionary is a
                          layer of code
- `chipping_rates::Array{Float,1}`: array of chipping rates in Hz for each code
- `channel::Complex{Int}`: `[NOT USED]` complex number used to define the
                           channel of the codes
    * set to `1+0im` if the channel is `I`
    * set to `0+1im` if the channel is `Q`
    * set to `1+1im` if the channel is `both`
- `include_codes::Array{Bool,1}`: array of `Bool` flags for each code to indicate
                                  whether a given code is being used
- `databits::Bool`: `Bool` flag to specify whether databits used in the set of codes
    * if `True`, databits are the last code in `CodeType.codes` array
- `similar_databits::Bool`: set to true if the databits are the same for all PRNs
"""
struct CodeType{T1,T2,T3,T4,T5,T6}
    name::String           # Name of signal type (i.e. "l1ca" or "l5i" or "l5q")
    code_num::T1           # Number of codes in signal type
                           # Includes primary, secondary, and nav codes
    codes::T2              # Vector of dictionaries
                           # Each dictionary is a primary, secondary, or
                           # or nav code which is specific for each PRN number
                           # The last dictionary in the vector of `codes` is
                           # assumed to be the navbits. This "code" will be
                           # ignored when generating replica signals for
                           # acquisition, fine acquisition, and signal tracking
                           # methods. Primary codes are considered the
                           # first code in `codes`.
    chipping_rates::T3     # Vector of chipping rates of length `code_num`
    code_lengths::T4       # Vector of code lengths of length `code_num`
    channel::T5            # Whether signal is on the I or Q channel
                           # If signal is on I, `channel = 1 + 0im`
                           # If signal is on Q, `channel = 0 + 1im`
                           # If signal is on both, `channel = 1 + 1im`
                           # `channel` is multiplied on resulting code value
                           # during signal generation
    include_codes::T6      # Bool vector of length `code_num`
                           # Used to determine if a given code is used for
                           # signal generation.
                           # Will be used in `calc_code_val` method to determine
                           # if a code value will be calculated and XOR`ed to
                           # other code values.
    databits::Bool         # Flag for if there are databits present in code
                           # Databits are assumed to be last entry in the
                           # vector of codes given to definecodetype funciton.
    similar_databits::Bool # True if databits are similar for all PRNs
end


"""
    SignalType

Holds common parms for a signal including its codes for each channel.

Fields:

- `name:String`: name of the signal type
- `I_codes::CodeType`: Custum codes for the I channel of the signal
    * set to `missing` if no `I_codes` given
- `Q_codes::CodeType`: Custum codes for the Q channel of the signal
    * set to `missing` if no `I_codes` given
- `sig_freq::Float64`: carrier frequency of the signal in Hz
- `B::Float64`: maximum bandwidth of all codes in Hz
- `include_I::Bool`: set to `True` if `I_codes` is given
- `include_Q::Bool`: set to `True` if `Q_codes` is given
"""
struct SignalType{T1,T2,T3,T4,T5}
    name::T1
    I_codes::T2
    Q_codes::T3
    sig_freq::T4
    B::T5
    include_I::Bool
    include_Q::Bool
end


"""
    definecodetype(codes::Vector, chipping_rates::Vector{Float64};
                   channel="both", databits=missing)

Defines a `CodeType` struct that represents the codes of a signal on a
given channel, `I`, `Q`, or `both`.

Arguments:

`codes::Vector`: array of codes, `not including databits`
    * the code is considered the primary code `(NOTE: this must be a dictionary)`
    * dictionary keys are the PRN numbers
    * all other codes after can be either dictionaries or arrays
    * if the code is an array, it is assumed that that layer of code is the
      same for all PRNs
`chipping_rates::Vector{Float}`: array of chipping rates for codes above, `not
                                 including databits`
`channel::String`: `[DEPRICATED]` the channel that the code will be on
    * can be either set to `I`, `Q`, or `both`
`databits::Vector (if given)`: if specified by user, must be a two element array
                               where the first element is the databits as either
                               a dictionary or array and the second element is
                               the chipping rate in Hz of the databits

Returns:

- `CodeType` structure
"""
function definecodetype(codes::Vector, chipping_rates::Vector{Float64};
                        channel="both", databits=missing)
    if ~isa(codes[1], Dict)
        error("Primary code must be a dictionary of codes. The primary code is the first index of `codes`.")
    end
    if ismissing(databits)
        code_num = length(codes)
        N = code_num
    else
        code_num = length(codes) + 1
        N = code_num - 1
    end
    if channel == "I"
        channel = 1 + 0im
        name = "I"
    elseif channel == "Q"
        channel = 0 + 1im
        name = "Q"
    elseif channel == "both"
        channel = 1 + 1im
        name = "IQ"
    else
        error("Invalid channel specified.")
    end
    code_keys = collect(keys(codes[1]))
    dict_type = Dict{eltype(code_keys),Vector{eltype(codes[1][code_keys[1]])}}
    signal_codes = Vector{dict_type}(undef, code_num)
    code_lengths = Array{Int}(undef, code_num)
    for i in 1:N
        if isa(codes[i], Dict)
            signal_codes[i] = codes[i]
            code_lengths[i] = length(codes[i][collect(keys(codes[i]))[1]])
        else
            code = dict_type()
            code_lengths[i] = length(codes[i])
            for code_key in code_keys
                code[code_key] = codes[i]
            end
            signal_codes[i] = code
        end
    end
    if ~ismissing(databits)
        push!(chipping_rates, databits[2])
        included_databits = true
        databit_codes = databits[1]
        if isa(databit_codes, Dict)
            similar_databits = false
            signal_codes[code_num] = databit_codes
            code_lengths[code_num] = length(databit_codes[collect(keys(databit_codes))[1]])
        else
            similar_databits = true
            code = dict_type()
            code_lengths[code_num] = length(databit_codes)
            for code_key in code_keys
                code[code_key] = databit_codes
            end
            signal_codes[code_num] = code
        end
    else
        included_databits = false
        similar_databits = false
    end
    include_codes = fill(true, code_num)
    return CodeType(name, code_num, signal_codes, chipping_rates, code_lengths,
                    channel, include_codes, included_databits, similar_databits)
end


"""
    definecodetype(code::Dict, chipping_rate::Float64,
                   channel="both"; databits=missing)

Defines a `CodeType` struct that represents the codes of a signal on a
given channel, `I`, `Q`, or `both`.

Arguments:

`code::Dict`: dictionary of primary codes where the keys are the PRN numbers
`chipping_rate::Vector{Float}`: array of chipping rates for codes above, `not
                                 including databits`
`channel::String`: `[DEPRICATED]` the channel that the code will be on
    * can be either set to `I`, `Q`, or `both`
`databits::Vector (if given)`: if specified by user, must be a two element array
                               where the first element is the databits as either
                               a dictionary or array and the second element is
                               the chipping rate in Hz of the databits

Returns:

- `CodeType` structure
"""
function definecodetype(code::Dict, chipping_rate::Float64,
                        channel="I"; databits=missing)
    if ~isa(code, Dict)
        error("Code given must be in the form of a dictionary.")
    end
    include_code = true
    if channel == "I"
        channel = 1 + 0im
        name = "I"
    elseif channel == "Q"
        name = "Q"
        channel = 0 + 1im
    elseif channel == "both"
        channel = 1 + 1im
        name = "IQ"
    else
        error("Invalid channel specified.")
    end
    code_keys = collect(keys(code))
    dict_type = Dict{eltype(code_keys),Vector{eltype(code[code_keys[1]])}}
    code_length = length(code[code_keys[1]])
    if ~ismissing(databits)
        code_num = 2
        signal_codes = Vector{dict_type}(undef, code_num)
        signal_codes[1] = code
        chipping_rate = [chipping_rate, databits[2]]
        included_databits = true
        databit_codes = databits[1]
        if isa(databit_codes, Dict)
            similar_databits = false
            signal_codes[code_num] = databit_codes
            code_length = [code_length, length(databit_codes[code_keys[1]])]
        else
            similar_databits = true
            code = dict_type()
            code_length = [code_length, length(databit_codes)]
            for code_key in code_keys
                code[code_key] = databit_codes
            end
            signal_codes[code_num] = code
        end
    else
        code_num = 1
        included_databits = false
        signal_codes = [code]
        similar_databits = false
    end
    include_codes = fill(true, code_num)
    return CodeType(name, code_num, signal_codes, chipping_rate, code_length,
                    channel, include_codes, included_databits, similar_databits)
end


"""
    definesignaltype(I_codes::CodeType, Q_codes::CodeType, sig_freq;
                     name="custom")

Defines a signal type with `CodeType` structs defined for both the `I` and `Q`
channels.

Arguments:

- `I_codes::CodeType`: Custum codes for the I channel of the signal
- `Q_codes::CodeType`: Custum codes for the Q channel of the signal
- `sig_freq::Float64`: carrier frequency of signal in Hz
- `name::String`: name of signal type (default is `custom`)
"""
function definesignaltype(I_codes::CodeType, Q_codes::CodeType, sig_freq;
                          name="custom")
    B_I = maximum(I_codes.chipping_rates)
    B_Q = maximum(Q_codes.chipping_rates)
    B = max(B_I, B_Q)
    return SignalType(name, I_codes, Q_codes, sig_freq, B, true, true)
end


"""
    definesignaltype(codes::CodeType, sig_freq, channel="I"; name="custom")

Defines a signal type with a `CodeType` struct defined for either the `I`, `Q`,
or `both` channels.

Arguments:

- `codes::CodeType`: Custum codes for either the `I`, `Q`, or `both` channels
- `sig_freq::Float64`: carrier frequency of signal in Hz
- `channel::String`: set to either `I`, `Q`, or `both` (default is `I`)
- `name::String`: name of signal type (default is `custom`)
"""
function definesignaltype(codes::CodeType, sig_freq, channel="I"; name="custom")
    B = maximum(codes.chipping_rates)
    if channel == "both"
        return SignalType(name, codes, codes, sig_freq, B, true, true)
    elseif channel == "I"
        return SignalType(name, codes, missing, sig_freq, B, true, false)
    elseif channel == "Q"
        return SignalType(name, missing, codes, sig_freq, B, false, true)
    else
        error("Invalid channel specified.")
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

Fields:

- `file_name::String`: name of data file
- `f_s::Float64`: sampling rate of data in Hz
- `f_if::Float64`: IF frequency of data in Hz
- `t_length::Float64`: length of data loaded in seconds
- `start_data_idx::Int64`: starting index in the data file
- `t::Array{Float64,1}`: time vector to accompany data
- `data::Array{Complex{Float64},1}`: loaded data from data file
- `data_type::String`: a `Symbol` type either `:sc4` or `:sc8` for 4 and 8 bit
                       comeplex data, respectively
- `data_start_time::T1`: `Tuple` of length 6
    * format is `(year, month, day, hour, minute, second)`
    * everything other than `second` must be an integer
    * `second` can be an integer or float
- `site_loc_lla::T2`: data collection site location
    * format is `(latitude, longitude, height)`
    * `latitude` and `longitude` are in degrees
    * `height` is in meters
- `sample_num::Int64`: number of data samples in `Data.data`
- `total_data_length::Float64`: total length of data file in seconds
- `nADC::Int64`: bit depth of data
    * `sc4` is 4 bit
    * `sc8` is 8 bit
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

A abstract struct for the replica signal structs. For use when custom signal
types are defined.

Fields:

- `name::String`
- `prn::Int64`
- `f_s::Float64`
- `t_length::Float64`
- `f_if::Float64`
- `f_d::Float64`
- `fd_rate::Float64`
- `Tsys::Float64`
- `CN0::Float64`
- `ϕ::Float64`
- `nADC::Int64`
- `code_start_idx::Float64`
- `init_code_phases_I::T1`
- `init_code_phases_Q::T2`
- `t::Array{Float64,1}`
- `data::Array{Complex{Float64},1}`
- `include_carrier::Bool`
- `include_adc::Bool`
- `include_thermal_noise::Bool`
- `include_databits_I::Bool`
- `include_phase_noise::Bool`
- `f_code_d_I::T2`
- `f_code_dd_I::T3`
- `f_code_d_Q::T4`
- `f_code_dd_Q::T5`
- `sample_num::Int64`
- `isreplica::Bool`
- `noexp::Bool`
- `thermal_noise::Array{Complex{Float64},1}`
- `phase_noise::Array{Float64,1}`
- `signal_type::T7`

"""
mutable struct ReplicaSignals{T1,T2,T3,T4,T5,T7} <: GNSSSignal
    name::String
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
    code_start_idx::Float64
    init_code_phases_I::T1
    init_code_phases_Q::T2
    t::Array{Float64,1}
    data::Array{Complex{Float64},1}
    include_carrier::Bool
    include_adc::Bool
    include_thermal_noise::Bool
    include_databits_I::Bool
    include_databits_Q::Bool
    include_phase_noise::Bool
    f_code_d_I::T2
    f_code_dd_I::T3
    f_code_d_Q::T4
    f_code_dd_Q::T5
    sample_num::Int64
    isreplica::Bool
    noexp::Bool
    thermal_noise::Array{Complex{Float64},1}
    phase_noise::Array{Float64,1}
    signal_type::T7
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
