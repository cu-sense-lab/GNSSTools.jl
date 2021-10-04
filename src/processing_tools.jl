"""
    process(signal::GNSSSignal, signal_type::SignalType, prn,
            channel="I"; σω=1000., fd_center=0., fd_range=5000.,
            cov_mult=1, q_a=1, q_mult=1, R_mult=1, dynamickf=true, 
            dll_b=5, state_num=3, fd_rate=0.,  figsize=missing, 
            saveto=missing, show_plot=true, fine_acq_method=:carrier,
            return_corrresult=false, use_fine_acq=true, σ_phi=π/2,
            h₀=1e-21, h₋₂=2e-20, acquisition_T=1e-3, fine_acq_T=10e-3, 
            tracking_T=1e-3, M=1, σᵩ²=missing)


Performs course acquisition and tracking on a `GNSSSignal` (either
`ReplicaSignal` or `GNSSData`).


Required Arguments:

- `signal::GNSSSignal`: simulated (`ReplicaSignal`) or real (`GNSSData`) data
- `signal_type`: `SignaType` struct that defines the signal being processed
    * i.e. L1 C/A or L5
- `prn`: PRN to processed (only one can be processed at a time)


Optional Arguments:

- `channel`: the channel to process (must be either `"I"` or `"Q"`)
             `(default = "I")`
- `σω`: Doppler rate uncertainty in Hz/s `(default = 1000)`
- `fd_center`: center frequency of the Doppler search area in Hz `(default = 0)`
- `fd_range`: range around the center frequency, `fd_center`, to search around
              during course acquisiton in Hz `(default = 5000)`
- `cov_mult`: scalar that can be applied to the covariance matrix in the
              tracking loop `(default = 1)`
- `q_a`: line of site platform dynamics in m²/s⁶ `(default = 1)`
- `q_mult`: scalar that can be applied to the process noise matrix, Q
            `(default = 1)`
- `dynamickf`: flag to specify whether dynamic KF gain is used or if only the
               steady state gain is used isntead `(default = true)`
- `dll_b`: delay lock loop bandwidth in Hz `(default = 5)`
- `state_num`: number of states to track in the carrier tracking loop
               `(default = 3)`
    * when `state_num = 2`: track only carrier phase and Doppler frequency
    * when `state_num = 3`: track the carrier phase, Doppler frequency, and
                            Doppler rate
- `fd_rate`: the initial estimate of the Doppler rate in Hz/s `(default = 0)`
- `figsize`: 2-element `Tuple` that specifies the height and width of the
             tracking results figure in inches `(default = missing)`
- `saveto`: path to save tracking results plot `(default = missing)`
- `show_plot`: flag to show tracking results plot `(default = true)`
- `fine_acq_method`: whether to use FFT or consecutive carrier based fine
                     acquisition method `(default = :carrier)`
    * can be either `:fft` of `:carrier`
- `return_corrresult`: flag, if set to true, will return the course correlation
                       result `(default = false)`
- `use_fine_acq`: flag to use a fine acquisition method `(default = true)`
- `acquisition_T`: integration time used for course acquisition in 
                  seconds `(default = 1e-3)`
- `fine_acq_T`: integration time used for fine acquisition in
                seconds `(default = 10e-3)`
- `tracking_T`: integration time used for code phase and carrier tracking in
                seconds `(default = 1e-3)`
- `M`: the number of coherent inetegrations to add non-coherently during
       course acquisition `(default = 1)`
- `σᵩ²`: the variance of the phase, also known as `R`
- `p_fa`: the probability of false alarm `(default = 1e-7)`
- `p_d`: the probability of detection `(default = 0.9)`


Returns:

- `results`: `FineAcquisitionResults` struct
- `trackresults`: `TrackResults` struct
- `corr_result`: `[OPTIONAL]` course correlation result returned if
                 `return_corrresult` is set to `true`
- `SNR_est`: `[OPTIONAL]` course correlation peak SNR returned if
             `return_corrresult` is set to `true`
"""
function process(signal::GNSSSignal, signal_type::SignalType, prn,
                 channel="I"; σω=1000., fd_center=0., fd_range=5000.,
                 cov_mult=1, q_a=1, q_mult=1, R_mult=1, dynamickf=true, 
                 dll_b=0.1, state_num=3, fd_rate=0.,  figsize=missing, 
                 saveto=missing, show_plot=true, fine_acq_method=:carrier,
                 return_corrresult=false, use_fine_acq=true, σ_phi=π/2,
                 h₀=1e-21, h₋₂=2e-20, acquisition_T=1e-3, fine_acq_T=10e-3, 
                 tracking_T=1e-3, M=1, σᵩ²=missing, return_Pd=false,
                 err_bin_num_f=0.25, fine_acq_sub_T=1e-3, p_fa=1e-7,
                 p_d=0.9)
    # Set up replica signals.
    f_s = signal.f_s
    if channel == "I"
        signal_type = definesignaltype(signal_type.I_codes,
                                       signal_type.sig_freq, "I")
    elseif channel == "Q"
        signal_type = definesignaltype(signal_type.Q_codes,
                                       signal_type.sig_freq, "Q")
    else
        error("Invalid channel `$(channel)`.")
    end
    replica = definereplica(signal_type, f_s, acquisition_T)
    # Perform course acquisition
    fd_course, n0_est, SNR_est, P_d, above_threshold,
    corr_result = courseacquisition(signal, 
                                    replica,
                                    prn; 
                                    fd_center=fd_center, 
                                    fd_range=fd_range,
                                    fd_rate=fd_rate,
                                    return_corrresult=true,
                                    M=M)
    # Estimate the C/N₀ and determine the value of the R matrix
    # since the estimate from the FFT based fine acquisition method
    # is always too large to use. This will only be used if the FFT
    # based fine acquisition method is used, or if no fine acquisition
    # method is used (when `use_fine_acq` = false)
    if ismissing(σᵩ²)
        processing_gain_est = processing_gain(acquisition_T, M; 
                                              p_fa=p_fa, p_d=p_d)
        CN0_est = SNR_est - processing_gain_est
        σᵩ² = phase_noise_variance(CN0_est, tracking_T)
    end
    # Perform fine acquisition using FFT based method
    # Returns structure containing the fine acquisition results,
    # fine Doppler (`fd_est`) and inital phase estimate (`phi_est`).
    # `results` structure also contains the initial uncertainties
    # for the states and measurements (`P` and `R`, respectively).
    # `σω` is the uncertainty of the Doppler rate. It is set to 1000Hz/s
    # by default.
    if use_fine_acq
        if fine_acq_method == :fft
            replicalong = definereplica(signal_type, f_s, fine_acq_T)
            results = fineacquisition(signal, replicalong, prn, fd_course,
                                      n0_est%replicalong.sample_num, 
                                      Val(fine_acq_method); σω=σω,
                                      fd_rate=fd_rate, 
                                      err_bin_num_f=err_bin_num_f)
        elseif fine_acq_method == :carrier
            M = floor(Int, fine_acq_T/fine_acq_sub_T)
            if M < 3
                error("Carrier based fine acquisition requires three iterations or more.")
            end
            replica = definereplica(signal_type, f_s, fine_acq_sub_T)
            results = fineacquisition(signal, replica, prn, fd_course,
                                      n0_est%replica.sample_num, 
                                      Val(fine_acq_method); σω=σω,
                                      fd_rate=fd_rate, M=M)              
        else
            error("Invalid value for argument `fine_acq_method`.")
        end
        # if fine_acq_method == :fft
        results.P[1,1] = σᵩ²
        # results.P[2,2] = (2π*err_bin_num_f/results.t_length)^2
        results.R[1] = σᵩ²
        # end
        R = results.R
        P = results.P
        phi_init = results.phi_init
        fd_est = results.fd_est
    else
        P = diagm([σᵩ², (2π*err_bin_num_f/acquisition_T)^2, (2π*σω)^2])
        R = [σᵩ²]
        phi_init = 0.
        fd_est = fd_course
        results = FineAcquisitionResults(prn, "N/A", fd_course, fd_rate, n0_est,
                                         acquisition_T, 0., fd_est, phi_init,
                                         "N/A", P, R)
    end
    # Peform tracking on signal using the initial estimates and
    # uncertainties calculated above.
    if tracking_T != acquisition_T 
        replica = definereplica(signal_type, f_s, tracking_T)
    end
    trackresults = trackprn(signal, replica, prn, phi_init,
                            fd_est, n0_est%replica.sample_num, 
                            P, R; DLL_B=dll_b,
                            state_num=state_num, dynamickf=dynamickf,
                            cov_mult=cov_mult, qₐ=q_a, q_mult=q_mult,
                            h₀=h₀, h₋₂=h₋₂, R_mult=R_mult, fd_rate=fd_rate)
    if show_plot
        plotresults(trackresults; saveto=saveto, figsize=figsize)
    end
    if return_corrresult && return_Pd
        return (results, trackresults, corr_result, SNR_est, 
                P_d, above_threshold)
    elseif return_corrresult && !return_Pd
        return (results, trackresults, corr_result, SNR_est)
    elseif !return_corrresult && return_Pd
        return (results, trackresults, P_d, above_threshold)
    else
        return (results, trackresults)
    end
end
