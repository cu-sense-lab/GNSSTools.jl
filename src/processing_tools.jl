"""
    process(signal::GNSSSignal, signal_type, prn, channel="both"; σω=1000.,
            fd_center=0., fd_range=5000., RLM=10, replica_t_length=1e-3,
            cov_mult=1, q_a=1, q_mult=1, dynamickf=true, dll_b=2,
            state_num=3, fd_rate=0., figsize=missing, saveto=missing,
            show_plot=true, fine_acq_method=:fft, M=10)


#


Required Arguments:

- `signal::GNSSSignal`:
- `signal_type`:
- `prn`:


Optional Arguments:

- `channel`:
- `σω`:
- `fd_center`:
- `fd_range`:
- `RLM`:
- `replica_t_length`:
- `cov_mult`:
- `q_a`:
- `q_mult`:
- `dynamickf`:
- `dll_b`:
- `state_num`:
- `fd_rate`:
- `figsize`:
- `saveto`:
- `show_plot`:
- `fine_acq_method`:
- `M`:
- `return_corrresult`: flag, if set to true, will return the course correlation
                       result


Returns:

- `results`:
- `trackresults`:
- `corr_result`: `[OPTIONAL]` course correlation result returned if
                 `return_corrresult` is set to `true`
"""
function process(signal::GNSSSignal, signal_type, prn, channel="both"; σω=1000.,
                 fd_center=0., fd_range=5000., RLM=10, replica_t_length=1e-3,
                 cov_mult=1, q_a=1, q_mult=1, dynamickf=true, dll_b=5,
                 state_num=3, fd_rate=0., figsize=missing, saveto=missing,
                 show_plot=true, fine_acq_method=:fft, M=10,
                 return_corrresult=false)
    # Set up replica signals. `replica_t_length` is used for
    # course acquisition and tracking, while `RLM*replica_t_length`
    # is used for fine acquisition only. The signal must be at least
    # as long as `RLM*replica_t_length`
    f_s = signal.f_s
    if isa(signal_type, SignalType)
        if channel == "both"
            # do nothing; keep as is
        elseif channel == "I"
            signal_type = definesignaltype(signal_type.I_codes,
                                           signal_type.sig_freq, "I")
        elseif channel == "Q"
            signal_type = definesignaltype(signal_type.Q_codes,
                                           signal_type.sig_freq, "Q")
        else
            error("Invalid channel `$(channel)`.")
        end
        replica = definesignal(signal_type, f_s, replica_t_length;
                               skip_noise_generation=true,
                               allocate_noise_vectors=false)
        replicalong = definesignal(signal_type, f_s, RLM*replica_t_length;
                                   skip_noise_generation=true,
                                   allocate_noise_vectors=false)
    else
        replica = definesignal(signal_type, f_s, replica_t_length)
        replicalong = definesignal(signal_type, f_s, RLM*replica_t_length)
    end
    if replicalong.sample_num > signal.sample_num
        error("Signal length equal to or greater than $(RLM*replica_t_length) seconds.")
    end
    # Initlalize 2D acquisition array
    Δfd = 1/replica.t_length  # Hz
    corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
    # Peform course acquisition
    courseacquisition!(corr_result, signal, replica, prn;
                       fd_center=fd_center, fd_range=fd_range,
                       fd_rate=fd_rate, Δfd=Δfd)
    # Calculate the code start index (`n0_est`), Doppler estimate (`fd_est`),
    # anc the SNR estimate (`SNR_est`)
    n0_est, fd_est, SNR_est = course_acq_est(corr_result, fd_center, fd_range,
                                             Δfd)
    # Perform fine acquisition using FFT based method
    # Returns structure containing the fine acquisition results,
    # fine Doppler (`fd_est`) and inital phase estimate (`phi_est`).
    # `results` structure also contains the initial uncertainties
    # for the states and measurements (`P` and `R`, respectively).
    # `σω` is the uncertainty of the Doppler rate. It is set to 1000Hz/s
    # by default.
    if fine_acq_method == :fft
        results = fineacquisition(signal, replicalong, prn, fd_est,
                                  n0_est, Val(fine_acq_method); σω=σω,
                                  fd_rate=fd_rate)
    elseif fine_acq_method == :carrier
        results = fineacquisition(signal, replica, prn, fd_est,
                                  n0_est, Val(fine_acq_method); σω=σω,
                                  fd_rate=fd_rate, M=RLM)
    else
        error("Invalid value for argument `fine_acq_method`.")
    end
    # Peform tracking on signal using the initial estimates and
    # uncertainties calculated above.
    trackresults = trackprn(signal, replica, prn, results.phi_init,
                            results.fd_est, results.n0_idx_course,
                            results.P, results.R; DLL_B=dll_b,
                            state_num=state_num, dynamickf=dynamickf,
                            cov_mult=cov_mult, qₐ=q_a, q_mult=q_mult)
    if show_plot
        plotresults(trackresults; saveto=saveto, figsize=figsize)
    end
    return (results, trackresults)
end
