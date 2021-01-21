##################################################################
######################## Generic Signal ##########################
##################################################################
"""
    definesignal(signal_type::SignalType, f_s, t_length; prn=1,
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., phi=0., nADC=4, include_carrier=true,
                 include_adc=true, include_thermal_noise=true,
                 code_start_idx=1., include_databits_I=true,
                 include_databits_Q=true, include_phase_noise=true,
                 phase_noise_scaler=1/10, name="custom",
                 skip_noise_generation=false, allocate_noise_vectors=true)


Define properties of locally generated generic signal
based off its type, PRN, etc.


Required Arguments:

- `signal_type::SignalType`: user defined signal type
    * see `definecodetype` and `definesignaltype`
- `f_s`: signal sampling frequency in Hz
- `t_length`: signal length in seconds


Optional Arguments:

- `prn`: satellite PRN number `(default = 1)`
- `f_if`: IF frequency in Hz `(default = 0Hz)`
- `f_d`: signal Doppler frequency in Hz `(default = 0Hz)`
- `fd_rate`: signal Doppler rate in Hz/s `(default = 0Hz/s)`
- `Tsys`: receiver noise temperature in Kelvin `(default = 535K)`
- `CN0`: signal carrier-to-noise ratio (C/N₀) `(default = 45dB⋅Hz)`
- `phi`: initial carrier phase in rad `(default = 0rad)`
- `nADC`: receiver bit depth `(default = 4 bits)`
- `include_carrier`: flag for modulating codes onto carrier `(default = true)`
- `include_adc`: flag for performing ADC quantization on signal `(default = true)`
- `include_thermal_noise`: flag for including thermal noise in signal generation
                           `(default = true)`
- `code_start_idx`: starting index in `ReplicaSignal.data` where all codes
                    start `(default = 1)`
- `include_databits_I`: flag for including I channel databits in simulated signal
                        `(default = true)`
- `include_databits_Q`: flag for including Q channel databits in simulated signal
                        `(default = true)`
- `include_phase_noise`: flag for including phase noise `(default = true)`
- `phase_noise_scaler`: standard deviation of phase noise in rad
                        `(default = 1/10rad)`
- `name`: name of signal `(default = custom)`
- `skip_noise_generation`: flag for skipping thermal and phase noise generation
                           `(default = false)`
- `allocate_noise_vectors`: flag for allocating memory for thermal and phase
                            noise vectors `(default = true)`
    * if set to `true`, thermal and phase noise vectors have the same length as
      `ReplicaSignal.data` vector
    * if set to `false`, thermal and phase noise vectors have length set to 0
      and the `include_thermal_noise` and `include_phase_noise` flags are set
      to `false`


Returns:

- `ReplicaSignal` struct
"""
function definesignal(signal_type::SignalType, f_s, t_length; prn=1,
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., phi=0., nADC=4, include_carrier=true,
                      include_adc=true, include_thermal_noise=true,
                      code_start_idx=1., include_databits_I=true,
                      include_databits_Q=true, include_phase_noise=true,
                      phase_noise_scaler=1/10, name="custom",
                      skip_noise_generation=false, allocate_noise_vectors=true)
    # Obtain sample number based off the sampling frequency and duration of
    # signal
    sample_num = Int(f_s * t_length)
    # Generate time vector
    t = calctvector(sample_num, f_s)
    # Calculate code chipping rates with Doppler applied for all codes on
    # each channel and their initial code phases
    sig_freq = signal_type.sig_freq
    I_codes = signal_type.I_codes
    Q_codes = signal_type.Q_codes
    # Get adjusted code chipping rates and chipping rate rates for I channel
    # codes
    if signal_type.include_I
        f_code_d_I = Array{Float64}(undef, I_codes.code_num)
        f_code_dd_I = Array{Float64}(undef, I_codes.code_num)
        init_code_phases_I = Array{Float64}(undef, I_codes.code_num)
        for i in 1:signal_type.I_codes.code_num
            # I channel
            f_code = I_codes.chipping_rates[i]
            code_length = I_codes.code_lengths[i]
            f_code_d, f_code_dd = calc_doppler_code_rate(f_code, sig_freq, f_d, fd_rate)
            f_code_d_I[i] = f_code_d
            f_code_dd_I[i] = f_code_dd
            init_code_phases_I[i] = calcinitcodephase(code_length, f_code_d,
                                                      f_code_dd, f_s, code_start_idx)
        end
        if I_codes.databits
            if ~include_databits_I
                I_codes.include_codes[end] = false
            else
                I_codes.include_codes[end] = true
            end
        end
    else
        # Set these vectors to zero length if there are no I channel codes
        f_code_d_I = Array{Float64}(undef, 0)
        f_code_dd_I = Array{Float64}(undef, 0)
        init_code_phases_I = Array{Float64}(undef, 0)
    end
    # Get adjusted code chipping rates and chipping rate rates for Q channel
    # codes
    if signal_type.include_Q
        f_code_d_Q = Array{Float64}(undef, Q_codes.code_num)
        f_code_dd_Q = Array{Float64}(undef, Q_codes.code_num)
        init_code_phases_Q = Array{Float64}(undef, Q_codes.code_num)
        for i in 1:signal_type.Q_codes.code_num
            # Q channel
            f_code = Q_codes.chipping_rates[i]
            code_length = Q_codes.code_lengths[i]
            f_code_d, f_code_dd = calc_doppler_code_rate(f_code, sig_freq, f_d, fd_rate)
            f_code_d_Q[i] = f_code_d
            f_code_dd_Q[i] = f_code_dd
            init_code_phases_Q[i] = calcinitcodephase(code_length, f_code_d,
                                                      f_code_dd, f_s, code_start_idx)
        end
        if Q_codes.databits
            if ~include_databits_Q
                Q_codes.include_codes[end] = false
            else
                Q_codes.include_codes[end] = true
            end
        end
    else
        # Set these vectors to zero length if there are no Q channel codes
        f_code_d_Q = Array{Float64}(undef, 0)
        f_code_dd_Q = Array{Float64}(undef, 0)
        init_code_phases_Q = Array{Float64}(undef, 0)
    end
    # Allocate space for signal. No simulated data is stored in yet. After
    # define the signal, user can use `generatesignal!(signal)` to generate
    # the I/Q samples which will be stored in `ReplicaSignal.data`.
    data = Array{Complex{Float64}}(undef, sample_num)
    isreplica = false
    noexp = false
    # Generate thermal noise and phase noise
    if allocate_noise_vectors
        if skip_noise_generation
            # Space for thermal and phase noise vectors are allocated. User can
            # use `definesignal!(signal; new_thermal_noise=true, new_phase_noise=true)`
            # to generate the noise in the vectors.
            thermal_noise = Array{Complex{Float64}}(undef, sample_num)
            phase_noise = Array{Float64}(undef, sample_num)
        else
            # Generate thermal and phase noise. Thermal noise is unscaled
            # but is scaled in `generatesignal!` when the signal is generated.
            # Phase noise is scaled based off the value of `phase_noise_scaler`.
            thermal_noise = randn(Complex{Float64}, sample_num)
            phase_noise = randn(Float64, sample_num)
            generate_phase_noise!(phase_noise, sample_num; scale=phase_noise_scaler)
        end
    else
        # Create 0 sized thermal and phase noise vectors and set the
        # `include_thermal_noise` and `include_phase_noise` flags to `false`.
        thermal_noise = Array{Complex{Float64}}(undef, 0)
        phase_noise = Array{Float64}(undef, 0)
        include_thermal_noise = false
        include_phase_noise = false
    end
    return ReplicaSignal(name, prn, f_s, t_length, f_if, f_d, fd_rate, Tsys,
                         CN0, phi, nADC, code_start_idx, init_code_phases_I,
                         init_code_phases_Q, t, data, include_carrier,
                         include_adc, include_thermal_noise, include_databits_I,
                         include_databits_Q, include_phase_noise, f_code_d_I,
                         f_code_dd_I, f_code_d_Q, f_code_dd_Q, sample_num,
                         isreplica, noexp, thermal_noise, phase_noise,
                         signal_type)
end


"""
    definesignal!(signal::ReplicaSignal;
                  prn=signal.prn, f_if=signal.f_if, f_d=signal.f_d,
                  fd_rate=signal.fd_rate, Tsys=signal.Tsys,
                  CN0=signal.CN0, phi=signal.phi, nADC=signal.nADC,
                  include_carrier=signal.include_carrier,
                  include_adc=signal.include_adc,
                  include_thermal_noise=signal.include_thermal_noise,
                  code_start_idx=signal.code_start_idx,
                  include_databits_I=signal.include_databits_I,
                  include_databits_Q=signal.include_databits_Q,
                  include_phase_noise=signal.include_phase_noise,
                  phase_noise_scaler=1/10, name=signal.name,
                  new_thermal_noise=false, new_phase_noise=false,
                  isreplica=signal.isreplica, noexp=signal.noexp,
                  new_databits=false)


Redefine properties of locally generated generic signal
based off its type, PRN, etc.


Required Arguments:

- `signal::ReplicaSignal`: the signal already defined by the user using
                            `definesignal`


Optional Arguments with Defaults Equal to `signal` Field Values:

- `prn`: satellite PRN number
- `f_if`: IF frequency in Hz
- `f_d`: signal Doppler frequency in Hz
- `fd_rate`: signal Doppler rate in Hz/s
- `Tsys`: receiver noise temperature in Kelvin
- `CN0`: signal carrier-to-noise ratio (C/N₀)
- `phi`: initial carrier phase in rad
- `nADC`: receiver bit depth
- `include_carrier`: flag for modulating codes onto carrier
- `include_adc`: flag for performing ADC quantization on signal
- `include_thermal_noise`: flag for including thermal noise in signal generation
- `code_start_idx`: starting index in `ReplicaSignal.data` where all codes
                    start
- `include_databits_I`: flag for including I channel databits in simulated signal
- `include_databits_Q`: flag for including Q channel databits in simulated signal
- `include_phase_noise`: flag for including phase noise
- `name`: name of signal
- `isreplica`: flag to identify whether signal is being for processing
    * replica signals being used for processing simulated or real signals
      should not include any noise source
    * this is usually set to `true` when being used by processing functions
      such as `courseacquisition!`, `fineacquisition`, and `trackprn`, only if
      this signal structure is being used as the replica signal, and is NOT
      the signal being processed
    * if `true`, a `ReplicaSignal` struct will not have noise added onto it
      and will not undergo ADC quantization
    * signal generation will be done using the second method of `generatesignal!`,
      `generatesignal!(signal::ReplicaSignal, isreplica::Bool)`
- `noexp`: used only if second method of `generatesignal!`, discussed
           above, is used
    * does not modulate codes onto carrier


Other Optional Arguments:

- `phase_noise_scaler`: standard deviation of phase noise in rad
                        `(default = 1/10rad)`
- `new_thermal_noise`: flag to generate new phase noise for signal
                       `(default = false)`
- `new_phase_noise`: flag to generate new phase noise for signal
                     `(default = false)`
- `new_databits`: flag to generate new databits for signal `(default = false)`


Modifies and Returns:

- `signal::ReplicaSignal`
"""
function definesignal!(signal::ReplicaSignal;
                       prn=signal.prn, f_if=signal.f_if, f_d=signal.f_d,
                       fd_rate=signal.fd_rate, Tsys=signal.Tsys,
                       CN0=signal.CN0, phi=signal.phi, nADC=signal.nADC,
                       include_carrier=signal.include_carrier,
                       include_adc=signal.include_adc,
                       include_thermal_noise=signal.include_thermal_noise,
                       code_start_idx=signal.code_start_idx,
                       include_databits_I=signal.include_databits_I,
                       include_databits_Q=signal.include_databits_Q,
                       include_phase_noise=signal.include_phase_noise,
                       phase_noise_scaler=1/10, name=signal.name,
                       new_thermal_noise=false, new_phase_noise=false,
                       isreplica=signal.isreplica, noexp=signal.noexp,
                       new_databits=false)
    # Calculate code chipping rates with Doppler applied for all codes on
    # each channel and their initial code phases
    f_s = signal.f_s
    signal_type = signal.signal_type
    sample_num = signal.sample_num
    sig_freq = signal_type.sig_freq
    I_codes = signal_type.I_codes
    Q_codes = signal_type.Q_codes
    # Get adjusted code chipping rates and chipping rate rates for I channel
    # codes
    if signal_type.include_I
        for i in 1:signal_type.I_codes.code_num
            # I channel
            f_code = I_codes.chipping_rates[i]
            code_length = I_codes.code_lengths[i]
            f_code_d, f_code_dd = calc_doppler_code_rate(f_code, sig_freq, f_d, fd_rate)
            signal.f_code_d_I[i] = f_code_d
            signal.f_code_dd_I[i] = f_code_dd
            signal.init_code_phases_I[i] = calcinitcodephase(code_length, f_code_d,
                                                             f_code_dd, f_s,
                                                             code_start_idx)
        end
        if I_codes.databits
            if ~include_databits_I
                signal.signal_type.I_codes.include_codes[end] = false
            else
                signal.signal_type.I_codes.include_codes[end] = true
            end
            # Generate new databits for I channel if `new_databits` flag is
            # `true`
            if new_databits
                prn_keys = collect(keys(I_codes.codes[1]))

                if I_codes.similar_databits
                    # If one array of databits was used for all PRNs, then the
                    # databit arrays for all PRNs are references to the same
                    # one. Therefore, changing one, will change them all.
                    rand!(signal.signal_type.I_codes.codes[end][prn_keys[1]], 0:1)
                else
                    # If a dictionary of databits was given instead, then each
                    # databit array is different. They are then separately
                    # regenerated.
                    for key in prn_keys
                        rand!(signal.signal_type.I_codes.codes[end][key], 0:1)
                    end
                end
            end
        end
    end
    # Get adjusted code chipping rates and chipping rate rates for Q channel
    # codes
    if signal_type.include_Q
        for i in 1:signal_type.Q_codes.code_num
            # Q channel
            f_code = Q_codes.chipping_rates[i]
            code_length = Q_codes.code_lengths[i]
            f_code_d, f_code_dd = calc_doppler_code_rate(f_code, sig_freq, f_d, fd_rate)
            signal.f_code_d_Q[i] = f_code_d
            signal.f_code_dd_Q[i] = f_code_dd
            signal.init_code_phases_Q[i] = calcinitcodephase(code_length, f_code_d,
                                                             f_code_dd, f_s,
                                                             code_start_idx)
        end
        # Generate new databits for Q channel if `new_databits` flag is
        # `true`
        if Q_codes.databits
            if ~include_databits_Q
                signal.signal_type.Q_codes.include_codes[end] = false
            else
                signal.signal_type.Q_codes.include_codes[end] = true
            end
            if new_databits
                prn_keys = collect(keys(Q_codes.codes[1]))
                if Q_codes.similar_databits
                    # If one array of databits was used for all PRNs, then the
                    # databit arrays for all PRNs are references to the same
                    # one. Therefore, changing one, will change them all.
                    rand!(signal.signal_type.Q_codes.codes[end][prn_keys[1]], 0:1)
                else
                    for key in prn_keys
                        # If a dictionary of databits was given instead, then each
                        # databit array is different. They are then separately
                        # regenerated.
                        rand!(signal.signal_type.Q_codes.codes[end][key], 0:1)
                    end
                end
            end
        end
    end
    # Generate new thermal noise and phase noise vectors if `new_thermal_noise`
    # and/or `new_phase_noise` flags are set to `true`
    if new_thermal_noise
        randn!(signal.thermal_noise)
    end
    if new_phase_noise
        randn!(signal.phase_noise)
        generate_phase_noise!(signal.phase_noise, sample_num; scale=phase_noise_scaler)
    end
    signal.name = name
    signal.prn = prn
    signal.f_if = f_if
    signal.f_d = f_d
    signal.fd_rate = fd_rate
    signal.Tsys = Tsys
    signal.CN0 = CN0
    signal.phi = phi
    signal.nADC = nADC
    signal.code_start_idx = code_start_idx
    signal.include_carrier = include_carrier
    signal.include_adc = include_adc
    signal.include_thermal_noise = include_thermal_noise
    signal.include_databits_I = include_databits_I
    signal.include_databits_Q = include_databits_Q
    signal.include_phase_noise = include_phase_noise
    signal.isreplica = isreplica
    signal.noexp = noexp
    return signal
end
