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
    descriminator::String
end


"""
    DLLParms

Struct for storing DLL parameters.
"""
struct DLLParms
    T::Float64
    B::Float64
    d::Int64
    descriminator::String
end


"""
    TrackResults

A struct containing the parameters for tracking and its results
"""
struct TrackResults{T1,T2}
    prn::Int64
    signaltype::String
    dll_parms::DLLParms
    pll_parms::PLLParms
    M::Int64
    T::Float64
    N::Int64
    data_file::String
    data_type::String
    data_nADC::Int64
    data_start_idx::Int64
    data_t_length::Float64
    data_total_t_length::Float64
    data_sample_num::Int64
    data_fs::Float64
    data_fif::Float64
    data_start_t::T1
    data_site_lla::T2
    data_init_n0::Float64
    data_init_code_chip::Float64
    data_init_phi::Float64
    data_init_fd::Float64
    t::Array{Float64,1}
    code_err_meas::Array{Float64,1}
    code_err_filt::Array{Float64,1}
    code_phase_meas::Array{Float64,1}
    code_phase_filt::Array{Float64,1}
    n0::Array{Float64,1}
    phi_measured::Array{Float64,1}
    phi_filtered::Array{Float64,1}
    phi::Array{Float64,1}
    delta_fd::Array{Float64,1}
    ZP::Array{Complex{Float64},1}
    SNR::Array{Float64,1}
    data_bits::Array{Int64,1}
    code_length::Int64
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
    b₂ = ωₙ^2*T^2 - 4*ξ*T*ωₙ
    descriminator = "ϕ = atan(Imag(ZP)/Real(ZP))"
    return PLLParms(T, ξ, B, ωₙ, a₀, a₁, a₂, b₀, b₁, b₂, descriminator)
end


"""
    definedll(T, B, d)

Define `DLLParms` struct.
"""
function definedll(T, B, d)
    descriminator = "Z4 = 1/d * (abs(ZE) - abs(ZL)) / (abs(ZE) + abs(ZL))"
    return DLLParms(T, B, d, descriminator)
end


"""
    measurephase(ZP)

Calculate the raw phase measurement.
"""
function measurephase(ZP)
    return atan((imag(ZP)/real(ZP)))
end


"""
    filtercarrierphase(pll_parms::PLLParms, ϕ_meas, ϕ_meas_1,
                       ϕ_meas_2, ϕ_filt_1, ϕ_filt_2)

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
    return 0.5 * (abs(ZE) - abs(ZL)) / (abs(ZE) + abs(ZL))
end


"""
    filtercodephase(dll_parms::DLLParms, current_code_err, last_filt_code_err)

Returns the filtered code phase measusurement.
"""
function filtercodephase(dll_parms::DLLParms, current_code_err, last_filt_code_err)
    return last_filt_code_err+4*dll_parms.T*dll_parms.B*(current_code_err-last_filt_code_err)
end


"""
    getcorrelatoroutput(data, replica, i, N, f_if, f_d, ϕ, d)

Calculate the early, prompt, and late correlator ouputs. Note that
replica already containts the prompt correlator. Be sure to set
the parameters to `replica` and run `generatesignal!(replica)` before
calling this method.
"""
function getcorrelatoroutput(data, replica, i, N, f_if, f_d, ϕ, d)
    ## Initialize correlator results
    # ze = 0. + 0im
    # zp = 0. + 0im
    # zl = 0. + 0im
    # datasegment = view(data.data, (i-1)*N+1:i*N)
    # ts = view(data.t, (i-1)*N+1:i*N)
    # ts = view(data.t, 1:N)
    ## Perform carrier and phase wipeoff and apply early, prompt, and late correlators
    # for j in 1:N
    #     @inbounds wipeoff = datasegment[j]*exp(-(2π*(f_if+f_d)*ts[j] + ϕ)*1im)
    #     @inbounds zp += conj(replica.data[j]) * wipeoff
    #     zeidx = j + d
    #     if zeidx > N
    #         zeidx = zeidx - N
    #     end
    #     zlidx = j - d
    #     if zlidx < 1
    #         zlidx = zlidx + N
    #     end
    #     @inbounds ze += conj(replica.data[zeidx]) * wipeoff
    #     @inbounds zl += conj(replica.data[zlidx]) * wipeoff
    # end
    wipeoff = data.data[(i-1)*N+1:i*N].*exp.(-(2π.*(f_if+f_d).*data.t[(i-1)*N+1:i*N] .+ ϕ).*1im)
    ze = sum(conj.(circshift(replica.data, -d)).*wipeoff)/N
    zp = sum(conj.(replica.data).*wipeoff)/N
    zl = sum(conj.(circshift(replica.data, d)).*wipeoff)/N
    # ze = ze/N
    # zp = zp/N
    # zl = zl/N
    return (ze, zp, zl)
end


"""
    trackprn(data, replica, prn, ϕ_init, fd_init, n0_init;
             DLL_B=5, PLL_B=15, damping=1.4, T=1e-3, M=1, d=2,
             t_length=data.t_length)

Perform code and phase tracking on data in `data`.

`replica` decides the signal type. Can pass optional arguments
that are minumum amount to track a given PRN.
"""
function trackprn(data, replica, prn, ϕ_init, fd_init, n0_idx_init;
                  DLL_B=5, PLL_B=15, damping=1.4, fd_rate=0.)
    # Assign signal specific parameters
    chipping_rate = replica.chipping_rate
    sig_freq = replica.sig_freq
    code_length = replica.code_length
    # Compute the spacing between the ZE, ZP, and ZL correlators
    d = Int64(floor(data.f_s/chipping_rate/2))
    # Initialize common variables and initial conditions
    T = replica.t_length
    t_length = data.t_length
    f_s = data.f_s
    f_if = data.f_if
    N = replica.sample_num
    f_d = fd_init
    ϕ = ϕ_init
    f_code_d = chipping_rate*(1. + f_d/sig_freq)
    n0_init = calcinitcodephase(code_length,
                                f_code_d, 0.,
                                f_s, n0_idx_init)
    n0 = n0_init
    M = Int64(floor(data.t_length/replica.t_length))
    t = Array(0:T:M*T-T)
    # Define DLL and PLL parameter structs
    dll_parms = definedll(T, DLL_B, d)
    pll_parms = definepll(T, PLL_B, damping)
    # Allocate array space for tracking results
    code_err_meas = Array{Float64}(undef, M)
    code_err_filt = Array{Float64}(undef, M)
    code_phase_meas = Array{Float64}(undef, M)
    code_phase_filt = Array{Float64}(undef, M)
    n0s = Array{Float64}(undef, M)
    phi_measured = Array{Float64}(undef, M)
    phi_filtered = Array{Float64}(undef, M)
    phi = Array{Float64}(undef, M)
    delta_fd = Array{Float64}(undef, M)
    ZP = Array{Complex{Float64}}(undef, M)
    SNR = Array{Float64}(undef, M)
    data_bits = Array{Int64}(undef, M)
    # Initialize 1ˢᵗ order DLL and 2ⁿᵈ PLL filters
    p = Progress(M, 1, "Tracking PRN $(prn)...")
    for i in 1:2
        # Calculate the current code start index
        t₀ = ((code_length-n0)%code_length)/f_code_d
        code_start_idx = t₀*f_s + 1
        # Set signal parameters
        definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, ϕ=0., f_if=0.,
                      include_carrier=true,
                      include_noise=false,
                      include_adc=false,
                      code_start_idx=code_start_idx,
                      isreplica=true, noexp=true)
        # Generate prompt correlator
        generatesignal!(replica)
        # Calculate early, prompt, and late correlator outputs
        ze, zp, zl = getcorrelatoroutput(data, replica, i, N, f_if, f_d, ϕ, d)
        # Estimate code phase error
        n0_err = Z4(dll_parms, ze, zp, zl)
        # Estimate carrier phase
        ϕ_meas = measurephase(zp)
        # Filter and update code phase
        if i > 1
            # Filter raw code phase error measurement
            n0_err_filtered = filtercodephase(dll_parms, n0_err, code_err_filt[i-1])
            # Calculate dfd
            dfd = (ϕ_meas - phi_measured[i-1])/(2π*T)
        else
            n0_err_filtered = n0_err
            dfd = 0.
        end
        # Save to allocated arrays
        code_err_meas[i] = n0_err
        code_err_filt[i] = n0_err_filtered
        code_phase_meas[i] = n0 + n0_err
        code_phase_filt[i] = n0 + n0_err_filtered
        n0s[i] = n0
        phi_measured[i] = ϕ_meas
        phi_filtered[i] = ϕ_meas
        phi[i] = ϕ
        delta_fd[i] = dfd
        ZP[i] = zp
        if real(zp) > 0
            data_bits[i] = 1
        else
            data_bits[i] = 0
        end
        # Calculate main code chipping rate at next `i`
        f_code_d = chipping_rate*(1. + f_d/sig_freq)
        # Update code phase with filtered code phase error and propagate to next `i`
        if i > 1
            n0 += n0_err_filtered# + f_code_d*T
        end
        # Updated and propagate carrier phase to next `i`
        ϕ += ϕ_meas# + (f_if + f_d)*T
        next!(p)
    end
    # Perform 1ˢᵗ and 2ⁿᵈ order DLL and PLL tracking, respectively
    for i in 3:M
        # Calculate the current code start index
        t₀ = ((code_length-n0)%code_length)/f_code_d
        code_start_idx = t₀*f_s + 1
        # Set signal parameters
        definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, ϕ=0., f_if=0.,
                      include_carrier=true,
                      include_noise=false,
                      include_adc=false,
                      code_start_idx=code_start_idx,
                      isreplica=true, noexp=true)
        # Generate prompt correlator
        generatesignal!(replica)
        # Calculate early, prompt, and late correlator outputs
        ze, zp, zl = getcorrelatoroutput(data, replica, i, N, f_if, f_d, ϕ, d)
        # Estimate code phase error
        n0_err = Z4(dll_parms, ze, zp, zl)
        # Estimate carrier phase
        ϕ_meas = measurephase(zp)
        # Filter carrier phase measurement
        ϕ_filt = filtercarrierphase(pll_parms, ϕ_meas,
                                    phi_measured[i-1], phi_measured[i-2],
                                    phi_filtered[i-1], phi_filtered[i-2])
        # Filter raw code phase error measurement
        n0_err_filtered = filtercodephase(dll_parms, n0_err, code_err_filt[i-1])
        # Calculate dfd
        dfd = (ϕ_meas - phi_measured[i-1])/(2π*T)
        # Save to allocated arrays
        code_err_meas[i] = n0_err
        code_err_filt[i] = n0_err_filtered
        code_phase_meas[i] = n0 + n0_err
        code_phase_filt[i] = n0 + n0_err_filtered
        n0s[i] = n0
        phi_measured[i] = ϕ_meas
        phi_filtered[i] = ϕ_filt
        phi[i] = ϕ
        delta_fd[i] = dfd
        ZP[i] = zp
        if real(zp) > 0
            data_bits[i] = 1
        else
            data_bits[i] = 0
        end
        # Calculate main code chipping rate at next `i`
        f_code_d = chipping_rate*(1. + f_d/sig_freq)
        # Update code phase with filtered code phase error and propagate to next `i`
        n0 += n0_err_filtered# + f_code_d*T
        # Update and propagate carrier phase to next `i`
        ϕ += ϕ_filt# + (f_if + f_d)*T
        next!(p)
    end
    # Return `TrackResults` struct
    # Specify what to store if `data` is a simulated signal struct
    # or `GNSSData` struct
    if typeof(data) == GNSSData
        file_name = data.file_name
        total_data_length = data.total_data_length
        data_type = data.data_type
        start_data_idx = data.start_data_idx
        site_loc_lla = data.site_loc_lla
        data_start_time = data.data_start_time
    else
        file_name = "N/A"
        total_data_length = data.t_length
        data_type = "N/A"
        start_data_idx = 1
        site_loc_lla = "N/A"
        data_start_time = "N/A"
    end
    return TrackResults(prn, gnsstypes[replica.type],
                        dll_parms, pll_parms, M, T, N,
                        file_name, data_type, data.nADC,
                        start_data_idx, data.t_length,
                        total_data_length, data.sample_num, f_s,
                        f_if, data_start_time, site_loc_lla,
                        float(n0_idx_init), n0_init, ϕ_init, fd_init, t,
                        code_err_meas, code_err_filt, code_phase_meas,
                        code_phase_filt, n0s, phi_measured, phi_filtered,
                        phi, delta_fd, ZP, SNR, data_bits, code_length)
end


"""
    plot(results::TrackResults, saveto=missing)
   
Plots the tracking results from the `trackprn` method. 
"""
function plotresults(results::TrackResults; saveto=missing)
    figure(figsize=(14,8))
    matplotlib.gridspec.GridSpec(3,2)
    # Plot code phase errors
    subplot2grid((3,2), (0,0), colspan=1, rowspan=1)
    plot(results.t, results.code_phase_meas.%results.code_length, "k.", label="Measured code phase")
    plot(results.t, results.code_phase_filt.%results.code_length, "b-", label="Filtered code phase")
    xlabel("Time (s)")
    ylabel("Code Phase (chips)")
    title("DLL Tracking")
    legend()
    # Plot filtered and measured phase errors
    subplot2grid((3,2), (0,1), colspan=1, rowspan=1)
    plot(results.t, results.phi_measured.*180 ./π, "k.", label="Measured ϕ")
    plot(results.t, results.phi_filtered.*180 ./π, "b-", label="Filtered ϕ")
    plot(results.t, results.phi_filtered.*180 ./π, "b.")
    xlabel("Time (s)")
    ylabel("ϕ (degrees)")
    # ylim([-180, 180])
    title("PLL Tracking")
    legend()
    subplot2grid((3,2), (1,0), colspan=2, rowspan=1)
    plot(results.t, results.delta_fd.+results.data_init_fd, "k.")
    xlabel("Time (s)")
    ylabel("Doppler (Hz)")
    title("Doppler Frequency Estimate")
    # subplot2grid((3,2), (1,1), colspan=1, rowspan=1)
    # plot(SNRs, "k.")
    # xlabel("Time (ms)")
    # ylabel("SNR (dB)")
    # title("Prompt Correlator SNR")
    # Plot ZP real and imaginary parts
    subplot2grid((3,2), (2,0), colspan=2, rowspan=1)
    plot(results.t, real(results.ZP), label="real(ZP)")
    plot(results.t, imag(results.ZP), label="imag(ZP)")
    xlabel("Time (s)")
    ylabel("ZP")
    title("Prompt Correlator Output")
    legend()
    suptitle("PRN $(results.prn)\nfd = $(Int(round(results.data_init_fd))) Hz\nn₀ = $(results.data_init_n0) samples")
    # subplots_adjust(hspace=0.4, wspace=0.4)
    subplots_adjust(hspace=0.4, wspace=0.2)
    if ~ismissing(saveto)
        savefig(saveto::String)
    end
end
