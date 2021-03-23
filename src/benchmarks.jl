"""
    CN0_monte_carlo(CN0, dopplers, doppler_rates, t_length, f_s, channel="I",
                    signal_type=define_l1ca_code_type(t_length); bins=1000,
                    iterations=100, prn=1, nADC=4, phase_noise_scaler=1/10,
                    include_phase_noise=true, include_thermal_noise=true,
                    include_databits_I=true, include_databits_Q=true,
                    f_if=0., include_adc=true, Tsys=535., include_carrier=true,
                    σω=1000., fd_range=5000., RLM=10, replica_t_length=1e-3,
                    cov_mult=1, q_a=1, q_mult=1, dynamickf=true, dll_b=5,
                    state_num=3, fine_acq_method=:fft, M=10)


Performs a benchmark test of signal simulation and processing by varying all
parameters, except the carrier-to-noise ratio, C/N₀, which is set by the user.


Required Arguments:

- `CN0`: carrier-to-noise ratio (C/N₀) of signal in dB
- `dopplers`: vector of Dopplers where duplicates are allowed
- `doppler_rates`: vector of Doppler rates where duplicates are allowed
- `t_length`: length of the signal to simulate and process in seconds
- `f_s`: sampling rate of signal in Hz


Optional Arguments:

- `channel`: channel to process signal on (must be either `"I"` or `"Q"`)
             `(default = "I")`
- `signal_type`: `SignalType` struct describing signal type
                 `(default = define_l1ca_code_type(t_length))`
- `bins`: number of bins that the Doppler vectors will be split into
          `(default = 1000)`
- `iterations`: number of times the simulation will be re-run and re-processed
                `(default = 100)`
- `prn`: PRN number to simulate `(default = 1)`
- `nADC`: bit depth of signal
- `phase_noise_scaler`: standard deviation of phase noise in rads
                       `(default = 1/10)`
- `include_phase_noise`: flag to include phase noise `(default = true)`
    * will be regenerated on each iteration
- `include_thermal_noise`: flag to include thermal noise `(default = true)`
    * will be regenerated on each iteration
- `include_databits_I`: flag to include databits on the I channel if they exist
                        `(default = true)`
- `include_databits_Q`: flag to include databits on the Q channel if they exist
                        `(default = true)`
- `f_if`: signal IF frequency in Hz `(default = 0Hz)`
- `include_adc`: flag to quantize signal `(default = true)`
- `Tsys`: system noise temperature in kelvin `(default = 535K)`
- `include_carrier`: flag to include carrier in signal `(default = true)`
- `σω`: initial uncertainty in Doppler rate `(default = 1000Hz/s)
- `fd_range`: range of Doppler above and below Doppler search center frequency
              in Hz `(default = 5000Hz)`
- `RLM`: a multiple of the replica signal length and used for FFT based fine
         acquisition `(default = 10)`
- `replica_t_length`: length of the replica signal in seconds `(default = 1e-3s)`
    * determines the length of the replica used in course acquisition and
      tracking
- `cov_mult`: scalar to inflate the initial covariance matrix `P₀`
              `(default = 1)`
- `q_a`: line of site platform dynamics in m²/s⁶ `(default = 1m²/s⁶)`
- `q_mult`: scalar to inflate the process noise matrix `Q` `(default = 1)`
- `dynamickf`: flag to specify if steady state KF gain is used or if KF can
               change from time step to time step `(default = true)`
- `dll_b`: bandwidth of the DLL filter in Hz `(default = 5Hz)`
- `state_num`: number of states to track `(default = 3)`
	* if set to `3`, carrier phase, Doppler, and Doppler rate are tracked
	* if set to `2`, only carrier phase and Doppler are tracked
- `fine_acq_method`:
- `M`: a multiple of the replica signal length and used for carrier based fine
       acquisition `(default = 10)`
- `h_parms`: 5 element vector containing oscillator h parameeters


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
                         cov_mult=1, q_a=1, q_mult=1, dynamickf=true, dll_b=5,
                         state_num=3, fine_acq_method=:fft, M=10, 
                         h_parms=h_parms_tcxo[1], σ_phi=π/2, fd_center=fd_center)
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
                          skip_noise_generation=true, receiver_h_parms=h_parms)
    ϕ = rand(0:0.0001:2π, iterations)
    code_start_idx = rand(1:signal.sample_num, iterations)
    T_num = floor(Int, signal.t_length/replica_t_length)
    t = calctvector(T_num, Δfd)
    p = Progress(iterations, 1, "Processing...")
    for i in 1:iterations
        f_d = sample(doppler_hist.edges[1], doppler_weights)
        if ismissing(fd_center)
            fd_center = round(f_d/Δfd)*Δfd
        end
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
                                           fd_center=fd_center,
                                           return_corrresult=false, 
                                           fine_acq=fine_acq, σ_phi=σ_phi)
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
