"""
    PLLParms

Struct for storing coefficients for the PLL filter. 
"""
struct PLLParms
    T::Float64
    damping::Float64
    B::Float64
    wn::Float64
    a0::Float64
    a1::Float64
    a2::Float64
    b0::Float64
    b1::Float64
    b2::Float64
end


"""
    DLLParms

Struct for storing DLL parameters.
"""
struct DLLParms
    T::Float64
    B::Float64
    d::Int64
end


"""
    TrackResults

A struct containing the parameters for tracking and its results
"""
struct TrackResults{T1,T2,T3,T4,T5,T6}
    prn::Int64
    signaltype::T1
    dll_parms::DLLParms
    pll_parms::PLLParms
    M::Int64
    integration_len::Float64
    integration_N::Int64
    data_file::String
    data_type::T2
    data_nADC::Int64
    data_start_idx::Int64
    data_t_length::Float64
    data_total_t_length::Float64
    data_sample_num::Int64
    data_fs::Float64
    data_if::Float64
    data_start_t::T3
    data_site_lla::T4
    data_init_n0::Float64
    data_init_code_chip::Float64
    data_init_phi::Float64
    data_init_fd::Float64
    t::Array{Float64,1}
    code_err_meas::Array{Float64,1}
    code_err_filt::Array{Float64,1}
    code_phase_meas::Array{Float64,1}
    code_phase_filtered::Array{Float64,1}
    phi_meassured::Array{Float64,1}
    phi_filtered::Array{Float64,1}
    delta_fd::Array{Float64,1}
    ZP::Array{Float64,1}
    SNR::Array{Float64,1}
    data_bits::Array{Int64,1}
    DLL_descriminator::T5
    PLL_descriminator::T6
end


"""
    definepll(T, B, damping)

Define `PLLParms` struct.
"""
function definepll(T, B, ξ)
    ωₙ = 4*B/(2*ξ + 1/(2*ξ))
    a₀ = ωₙ^2*T^2 + 4*ξ*T*ωₙ + 4
    a₁ = 2*ωₙ^2*T^2 - 8
    a₂ = ωₙ^2*T^2 - 4*ξ*T*ωₙ + 4
    b₀ = ωₙ^2*T^2 + 4*ξ*T*ωₙ
    b₁ = 2*ωₙ^2*T^2
    b₂ = ωₙ^2*T^2 - 4*ξ*ωₙ*T
    return PLLParms(T, ξ, B, ωₙ, a₀, a₁, a₂, b₀, b₁, b₂)
end


"""
    definedll(T, B, d)

Define `DLLParms` struct.
"""
function definedll(T, B, d)
    return DLLParms(T, B, d)
end


"""
    meassurephase(ZP)

Calculate the raw phase meassurement.
"""
function meassurephase(ZP)
    return atan((imag(ZP)/real(ZP)))
end


"""
    filtercarrierphase(pll_parms::PLLParms, ϕ_meas,
                       ϕ_filt_1, ϕ_filt_2)

Fiters the raw phase measurement using 2ⁿᵈ order PLL filter.
"""
function filtercarrierphase(pll_parms::PLLParms, ϕ_meas, ϕ_meas_1,
                            ϕ_meas_2, ϕ_filt_1, ϕ_filt_2)
    a₀ = pll_parms.a0
    a₁ = pll_parms.a1
    a₂ = pll_parms.a2
    b₀ = pll_parms.b0
    b₁ = pll_parms.b1
    b₂ = pll_parms.b2
    return (b₀*ϕ_meas + b₁*ϕ_meas_1 + b₂*ϕ_meas_2 - a₁*ϕ_filt_1 - a₂*ϕ_filt_2)/a₀
end


"""
    Z4(dll_parms::DLLParms, ZE, ZP, ZL)

Calculates the code phase error.
"""
function Z4(dll_parms::DLLParms, ZE, ZP, ZL)
    return 1/dll_parms.d * (abs(ZE) - abs(ZL)) / (abs(ZE) + abs(ZL))
end


"""
    filtercodephase(dll_parms::DLLParms, current_code_err, last_code_err)

Returns the filtered code phase measusurement.
"""
function filtercodephase(dll_parms::DLLParms, current_code_err, last_filt_code_err)
    return last_filt_code_err+4*dll_parms.T*dll_parms*B*(current_code_err-last_filt_code_err)
end


"""
    trackprn(data::GNSSData, replica, results=missing)

Perform code and phase tracking on data in `data`.

`replica` decides the signal type and `results` decides
the tracking parameters. Can also pass optional arguments
that are minumum amount to constuct `TrackResults` struct. 
"""
function trackprn(data::GNSSData, replica, prn, ϕ_init, fd_init, n0_init;
                  DLL_B=5, PLL_B=15, damping=1.4, integration_len=1e-3, M=1)
    
end