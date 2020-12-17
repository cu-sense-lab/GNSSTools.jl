"""
    CN0_monte_carlo(CN0, dopplers, doppler_rates, t_length, f_s, channel="I",
                    signal_type=define_l1ca_code_type(t_length); bins=1000,
                    iterations=100, prn=1, nADC=4, phase_noise_scaler=1/10,
                    include_phase_noise=true, include_thermal_noise=true,
                    include_databits_I=true, include_databits_Q=true,
                    f_if=0., include_adc=true, Tsys=535., include_carrier=true,
                    σω=1000., fd_range=5000., RLM=10, replica_t_length=1e-3,
                    cov_mult=1, q_a=1, q_mult=1, dynamickf=true, dll_b=2,
                    state_num=3, fine_acq_method=:fft, M=10)


#


Required Arguments:

- `CN0`:
- `dopplers`:
- `doppler_rates`:
- `t_length`:
- `f_s`:


Optional Arguments:

- `channel`:
- `signal_type`:
- `bins`:
- `iterations`:
- `prn`:
- `nADC`:
- `phase_noise_scaler`:
- `include_phase_noise`:
- `include_thermal_noise`:
- `include_databits_I`:
- `include_databits_Q`:
- `f_if`:
- `include_adc`:
- `Tsys`:
- `include_carrier`:
- `σω`:
- `fd_range`:
- `RLM`:
- `replica_t_length`:
- `cov_mult`:
- `q_a`:
- `q_mult`:
- `dynamickf`:
- `dll_b`:
- `state_num`:
- `fine_acq_method`:
- `M`:


Returns:

-
"""
function CN0_monte_carlo(CN0, dopplers, doppler_rates, t_length, f_s, channel="I",
                         signal_type=define_l1ca_code_type(t_length); bins=1000,
                         iterations=100, prn=1, nADC=4, phase_noise_scaler=1/10,
                         include_phase_noise=true, include_thermal_noise=true,
                         include_databits_I=true, include_databits_Q=true,
                         f_if=0., include_adc=true, Tsys=535., include_carrier=true,
                         σω=1000., fd_range=5000., RLM=10, replica_t_length=1e-3,
                         cov_mult=1, q_a=1, q_mult=1, dynamickf=true, dll_b=2,
                         state_num=3, fine_acq_method=:fft, M=10)
    if (channel != "I") && (channel != "Q")
        error("Invalid channel. Channel to be processed must be either I or Q.")
    end
    doppler_hist = fit(Histogram, dopplers, nbins=bins)
    doppler_weights = Weights(doppler_hist.weights)
    doppler_rate_hist = fit(Histogram, doppler_rates, nbins=bins)
    doppler_rate_weights = Weights(doppler_rate_hist.weights)
    Δfd = 1/replica_t_length
    if signal_type.include_I
        include_databits_I = signal_type.I_codes.databits
    else
        include_databits_I = false
    end
    if signal_type.include_Q
        include_databits_Q = signal_type.Q_codes.databits
    else
        include_databits_Q = false
    end
    include_databits = include_databits_I || include_databits_Q
    signal = definesignal(signal_type, f_s, t_length; prn=prn, CN0=CN0,
                          nADC=nADC, f_if=f_if, include_adc=include_adc,
                          Tsys=Tsys, include_carrier=include_carrier,
                          phase_noise_scaler=phase_noise_scaler,
                          include_phase_noise=include_phase_noise,
                          include_thermal_noise=include_thermal_noise,
                          include_databits_I=include_databits_I,
                          include_databits_Q=include_databits_Q,
                          skip_noise_generation=true)
    ϕ = rand(0:0.0001:2π, iterations)
    code_start_idx = rand(1:signal.sample_num, iterations)
    T_num = floor(Int, signal.t_length/replica_t_length)
    t = calctvector(T_num, Δfd)
    p = Progress(iterations, 1, "Processing...")
    for i in 1:iterations
        f_d = sample(doppler_hist.edges[1], doppler_weights)
        fd_center = fd_center = round(f_d/Δfd)*Δfd
        fd_rate = sample(doppler_rate_hist.edges[1], doppler_rate_weights)
        definesignal!(signal; f_d=f_d, fd_rate=fd_rate, phi=ϕ[i],
                      code_start_idx=code_start_idx[i], new_databits=true,
                      new_thermal_noise=true, new_phase_noise=true,)
        generatesignal!(signal)
        acqresults, trackresults = process(signal, signal_type, prn, channel;
                                           show_plot=false, σω=σω, q_a=q_a,
                                           fd_range=fd_range, RLM=RLM,
                                           q_mult=q_mult, cov_mult=cov_mult,
                                           dynamickf=dynamickf, dll_b=dll_b,
                                           state_num=state_num, fd_rate=fd_rate,
                                           fine_acq_method=fine_acq_method, M=M,
                                           replica_t_length=replica_t_length,
                                           fd_center=fd_center)
        # Truth data
        if channel == "I"
            init_code_phase = calcinitcodephase(signal_type.I_codes.code_lengths[1],
                                                signal.f_code_d_I[1],
                                                signal.f_code_dd_I[1],
                                                f_s, code_start_idx[i])
            truth_code = calccodeidx.(init_code_phase, signal.f_code_d_I[1],
                                      signal.f_code_dd_I[1], t,
                                      signal_type.I_codes.code_lengths[1])
            if signal.signal_type.I_codes.databits
                truth_databits = signal.signal_type.I_codes.codes[end][prn]
            else
                truth_databits = missing
            end
        elseif channel == "Q"
            init_code_phase = calcinitcodephase(signal_type.Q_codes.code_lengths[1],
                                                signal.f_code_d_Q[1],
                                                signal.f_code_dd_Q[1],
                                                f_s, code_start_idx[i])
            truth_code = calccodeidx.(init_code_phase, signal.f_code_d_Q[1],
                                      signal.f_code_dd_Q[1], t,
                                      signal_type.Q_codes.code_lengths[1])
            if signal.signal_type.Q_codes.databits
                truth_databits = signal.signal_type.Q_codes.codes[end][prn]
            else
                truth_databits = missing
            end
        else
            error("Invalid channel specified. Channel must be either I or Q.")
        end
        truth_doppler = f_d .+ fd_rate.*t
        next!(p)
    end
end
