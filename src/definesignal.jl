##################################################################
######################## Generic Signal ##########################
##################################################################
"""
    definesignal(signal_type::SignalType, f_s, t_length; prn=1,
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., ϕ=0., nADC=4, include_carrier=true,
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
- `code_start_idx`: starting index in `ReplicaSignals.data` where all codes
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
      `ReplicaSignals.data` vector
    * if set to `false`, thermal and phase noise vectors have length set to 0
      and the `include_thermal_noise` and `include_phase_noise` flags are set
      to `false`


Returns:

- `ReplicaSignals` struct
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
    # Allocate space for signal
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
        include_phase_noise = true
    end
    return ReplicaSignals(name, prn, f_s, t_length, f_if, f_d, fd_rate, Tsys,
                          CN0, phi, nADC, code_start_idx, init_code_phases_I,
                          init_code_phases_Q, t, data, include_carrier,
                          include_adc, include_thermal_noise, include_databits_I,
                          include_databits_Q, include_phase_noise, f_code_d_I,
                          f_code_dd_I, f_code_d_Q, f_code_dd_Q, sample_num,
                          isreplica, noexp, thermal_noise, phase_noise,
                          signal_type)
end


"""
    definesignal!(signal::ReplicaSignals;
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

- `signal::ReplicaSignals`: the signal already defined by the user using
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
- `code_start_idx`: starting index in `ReplicaSignals.data` where all codes
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
    * if `true`, a `ReplicaSignals` struct will not have noise added onto it
      and will not undergo ADC quantization
    * signal generation will be done using the second method of `generatesignal!`,
      `generatesignal!(signal::ReplicaSignals, isreplica::Bool)`
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

- `signal::ReplicaSignals`
"""
function definesignal!(signal::ReplicaSignals;
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
            if new_databits
                prn_keys = collect(keys(I_codes.codes[1]))
                if I_codes.similar_databits
                    rand!(signal.signal_type.I_codes.codes[end][prn_keys[1]], 0:1)
                else
                    for key in prn_keys
                        rand!(signal.signal_type.I_codes.codes[end][key], 0:1)
                    end
                end
            end
        end
    end
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
        if Q_codes.databits
            if ~include_databits_Q
                signal.signal_type.Q_codes.include_codes[end] = false
            else
                signal.signal_type.Q_codes.include_codes[end] = true
            end
            if new_databits
                prn_keys = collect(keys(Q_codes.codes[1]))
                if Q_codes.similar_databits
                    rand!(signal.signal_type.Q_codes.codes[end][prn_keys[1]], 0:1)
                else
                    for key in prn_keys
                        rand!(signal.signal_type.Q_codes.codes[end][key], 0:1)
                    end
                end
            end
        end
    end
    # Generate thermal noise and phase noise
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


#------------------------------------------------------------------------------
#                             OLD IMPLEMENTATON
#------------------------------------------------------------------------------


##################################################################
######################### L1 C/A Signal ##########################
##################################################################
"""
    definesignal(type::Val{:l1ca}, f_s, t_length; prn=1,
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., ϕ=0., nADC=4, B=2.046e6,
                 include_carrier=true, include_adc=true,
                 include_noise=true, code_start_idx=1,
                 include_databits=true, sig_freq=missing)

Define properties of locally generated L1 C/A signal
based off its type, PRN, etc.
"""
function definesignal(type::Val{:l1ca}, f_s, t_length; prn=1,
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=4, B=2.046e6,
                      include_carrier=true, include_adc=true,
                      include_noise=true, code_start_idx=1,
                      include_databits=true, sig_freq=missing)
    if ismissing(sig_freq)
        sig_freq = L1_freq
    end
    sample_num = Int(f_s * t_length)
    # Generate time vector
    t = calctvector(sample_num, f_s)
    # Calculate number of data bits (use `ceil` to ensure that databits do not repeat)
    db_code_length = Int64(ceil(l1ca_db_chipping_rate*t_length))
    # Calculate code chipping rates with Doppler applied
    # L1 C/A
    f_l1ca_d = l1ca_chipping_rate*(1. + f_d/sig_freq)
    f_l1ca_dd = l1ca_chipping_rate*fd_rate/sig_freq
    # Data bit sequence
    f_db_d = l1ca_db_chipping_rate*(1. + f_d/sig_freq)
    f_db_dd = l1ca_db_chipping_rate*fd_rate/sig_freq
    # Calculate the L1 C/A and data bit code phase offsets
    l1ca_init_code_phase = calcinitcodephase(l1ca_code_length,
                                             f_l1ca_d, f_l1ca_dd,
                                             f_s, code_start_idx)
    db_init_code_phase = calcinitcodephase(db_code_length,
                                           f_db_d, f_db_dd,
                                           f_s, code_start_idx)
    # Allocate space for signal
    data = Array{Complex{Float64}}(undef, sample_num)
    isreplica = false
    noexp = false
    # Generate random databit vector
    databits = rand(0:1, db_code_length)
    return L1CASignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                      Tsys, CN0, ϕ, nADC, B, code_start_idx,
                      l1ca_init_code_phase, db_init_code_phase, t,
                      data, include_carrier, include_adc,
                      include_noise, include_databits,
                      f_l1ca_d, f_l1ca_dd,
                      f_db_d, f_db_dd, sample_num,
                      isreplica, noexp, l1ca_chipping_rate,
                      sig_freq, l1ca_code_length, db_code_length, databits)
end


"""
    definesignal!(signal::L1CASignal;
                  prn=signal.prn, f_d=signal.f_d,
                  f_if=signal.f_if, fd_rate=signal.fd_rate,
                  Tsys=signal.Tsys, CN0=signal.CN0,
                  ϕ=signal.ϕ, nADC=signal.nADC,
                  B=signal.B,
                  include_carrier=signal.include_carrier,
                  include_adc=signal.include_adc,
                  include_noise=signal.include_noise,
                  code_start_idx=signal.code_start_idx,
                  isreplica=signal.isreplica,
                  noexp=signal.noexp,
                  include_databits=signal.include_databits,
                  sig_freq=signal.sig_freq)

Redefine properties of locally generated L5Q signal
based off its type, PRN, etc.

Sampling rate (f_s) and signal length (t_length)
are the only parameters that cannot be redefined.

"""
function definesignal!(signal::L1CASignal;
                       prn=signal.prn, f_d=signal.f_d,
                       f_if=signal.f_if, fd_rate=signal.fd_rate,
                       Tsys=signal.Tsys, CN0=signal.CN0,
                       ϕ=signal.ϕ, nADC=signal.nADC,
                       B=signal.B,
                       include_carrier=signal.include_carrier,
                       include_adc=signal.include_adc,
                       include_noise=signal.include_noise,
                       code_start_idx=signal.code_start_idx,
                       isreplica=signal.isreplica,
                       noexp=signal.noexp,
                       include_databits=signal.include_databits,
                       sig_freq=signal.sig_freq)
    ## Calculate code chipping rates with Doppler applied
    # L1 C/A
    f_l1ca_d = l1ca_chipping_rate*(1. + f_d/sig_freq)
    f_l1ca_dd = l1ca_chipping_rate*fd_rate/sig_freq
    # Data bit sequence
    f_db_d = l1ca_db_chipping_rate*(1. + f_d/sig_freq)
    f_db_dd = l1ca_db_chipping_rate*fd_rate/sig_freq
    # Calculate the L1 C/A and data bit code phase offsets
    l1ca_init_code_phase = calcinitcodephase(l1ca_code_length,
                                            f_l1ca_d, f_l1ca_dd,
                                            signal.f_s,
                                            code_start_idx)
    db_init_code_phase = calcinitcodephase(signal.db_code_length,
                                           f_db_d, f_db_dd,
                                           signal.f_s,
                                           code_start_idx)
    # Store udated variables to `L1CASignal` struct
    signal.prn = prn
    signal.f_d = f_d
    signal.f_if = f_if
    signal.fd_rate = fd_rate
    signal.Tsys = Tsys
    signal.CN0 = CN0
    signal.ϕ = ϕ
    signal.nADC = nADC
    signal.B = B
    signal.include_carrier = include_carrier
    signal.include_adc = include_adc
    signal.include_noise = include_noise
    signal.f_l1ca_d = f_l1ca_d
    signal.f_l1ca_dd = f_l1ca_dd
    signal.f_db_d = f_db_d
    signal.f_db_dd = f_db_dd
    signal.l1ca_init_code_phase = l1ca_init_code_phase
    signal.db_init_code_phase = db_init_code_phase
    signal.code_start_idx = code_start_idx
    signal.isreplica = isreplica
    signal.noexp = noexp
    signal.include_databits = include_databits
    signal.sig_freq = sig_freq
    return signal
end


##################################################################
########################### L5Q Signal ###########################
##################################################################
"""
    definesignal(type::Val{:l5q}, f_s, t_length; prn=1,
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., ϕ=0., nADC=4, B=2.046e7,
                 include_carrier=true, include_adc=true,
                 include_noise=true, code_start_idx=1,
                 sig_freq=missing)

Define properties of locally generated L5Q signal
based off its type, PRN, etc.
"""
function definesignal(type::Val{:l5q}, f_s, t_length; prn=1,
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=4, B=2.046e7,
                      include_carrier=true, include_adc=true,
                      include_noise=true, code_start_idx=1,
                      sig_freq=missing)
    if ismissing(sig_freq)
        sig_freq = L5_freq
    end
    sample_num = Int(f_s * t_length)
    # Generate time vector
    t = calctvector(sample_num, f_s)
    # Calculate code chipping rates with Doppler applied
    # L5Q
    f_l5q_d = L5_chipping_rate*(1. + f_d/sig_freq)
    f_l5q_dd = L5_chipping_rate*fd_rate/sig_freq
    # Neuman sequence
    f_nh_d = nh_chipping_rate*(1. + f_d/sig_freq)
    f_nh_dd = nh_chipping_rate*fd_rate/sig_freq
    # Calculate the L5Q and nh code phase offsets
    l5q_init_code_phase = calcinitcodephase(L5_code_length,
                                            f_l5q_d, f_l5q_dd,
                                            f_s, code_start_idx)
    nh_init_code_phase = calcinitcodephase(nh20_code_length,
                                           f_nh_d, f_nh_dd,
                                           f_s, code_start_idx)
    # Allocate space for signal
    data = Array{Complex{Float64}}(undef, sample_num)
    isreplica = false
    noexp = false
    return L5QSignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                     Tsys, CN0, ϕ, nADC, B, code_start_idx,
                     l5q_init_code_phase, nh_init_code_phase, t,
                     data, include_carrier, include_adc,
                     include_noise,
                     f_l5q_d, f_l5q_dd,
                     f_nh_d, f_nh_dd, sample_num,
                     isreplica, noexp, L5_chipping_rate,
                     sig_freq, L5_code_length)
end


"""
    definesignal!(signal::L5QSignal;
                  prn=signal.prn, f_d=signal.f_d,
                  f_if=signal.f_if, fd_rate=signal.fd_rate,
                  Tsys=signal.Tsys, CN0=signal.CN0,
                  ϕ=signal.ϕ, nADC=signal.nADC,
                  B=signal.B,
                  include_carrier=signal.include_carrier,
                  include_adc=signal.include_adc,
                  include_noise=signal.include_noise,
                  code_start_idx=signal.code_start_idx,
                  isreplica=signal.isreplica,
                  noexp=signal.noexp,
                  sig_freq=signal.sig_freq)

Redefine properties of locally generated L5Q signal
based off its type, PRN, etc.

Sampling rate (f_s) and signal length (t_length)
are the only parameters that cannot be redefined.

"""
function definesignal!(signal::L5QSignal;
                       prn=signal.prn, f_d=signal.f_d,
                       f_if=signal.f_if, fd_rate=signal.fd_rate,
                       Tsys=signal.Tsys, CN0=signal.CN0,
                       ϕ=signal.ϕ, nADC=signal.nADC,
                       B=signal.B,
                       include_carrier=signal.include_carrier,
                       include_adc=signal.include_adc,
                       include_noise=signal.include_noise,
                       code_start_idx=signal.code_start_idx,
                       isreplica=signal.isreplica,
                       noexp=signal.noexp,
                       sig_freq=signal.sig_freq)
    ## Calculate code chipping rates with Doppler applied
    # L5Q
    f_l5q_d = L5_chipping_rate*(1. + f_d/sig_freq)
    f_l5q_dd = L5_chipping_rate*fd_rate/sig_freq
    # Neuman sequence
    f_nh_d = nh_chipping_rate*(1. + f_d/sig_freq)
    f_nh_dd = nh_chipping_rate*fd_rate/sig_freq
    # Calculate the L5Q and nh code phase offsets
    l5q_init_code_phase = calcinitcodephase(L5_code_length,
                                            f_l5q_d, f_l5q_dd,
                                            signal.f_s,
                                            code_start_idx)
    nh_init_code_phase = calcinitcodephase(nh20_code_length,
                                           f_nh_d, f_nh_dd,
                                           signal.f_s,
                                           code_start_idx)
    # Store udated variables to "L5QSignal" struct
    signal.prn = prn
    signal.f_d = f_d
    signal.f_if = f_if
    signal.fd_rate = fd_rate
    signal.Tsys = Tsys
    signal.CN0 = CN0
    signal.ϕ = ϕ
    signal.nADC = nADC
    signal.B = B
    signal.include_carrier = include_carrier
    signal.include_adc = include_adc
    signal.include_noise = include_noise
    signal.f_l5q_d = f_l5q_d
    signal.f_l5q_dd = f_l5q_dd
    signal.f_nh_d = f_nh_d
    signal.f_nh_dd = f_nh_dd
    signal.l5q_init_code_phase = l5q_init_code_phase
    signal.nh_init_code_phase = nh_init_code_phase
    signal.code_start_idx = code_start_idx
    signal.isreplica = isreplica
    signal.noexp = noexp
    signal.sig_freq = sig_freq
    return signal
end


##################################################################
########################### L5I Signal ###########################
##################################################################
"""
    definesignal(type::Val{:l5i}, f_s, t_length; prn=1,
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., ϕ=0., nADC=4, B=2.046e7,
                 include_carrier=true, include_adc=true,
                 include_noise=true, code_start_idx=1,
                 include_databits=true, sig_freq=missing)

Define properties of locally generated L5I signal
based off its type, PRN, etc.
"""
function definesignal(type::Val{:l5i}, f_s, t_length; prn=1,
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=4, B=2.046e7,
                      include_carrier=true, include_adc=true,
                      include_noise=true, code_start_idx=1,
                      include_databits=true, sig_freq=missing)
    if ismissing(sig_freq)
        sig_freq = L5_freq
    end
    sample_num = Int(f_s * t_length)
    # Generate time vector
    t = calctvector(sample_num, f_s)
    # Calculate number of data bits (use `ceil` to ensure that databits do not repeat)
    db_code_length = Int64(ceil(L5_db_chipping_rate*t_length))
    # Calculate code chipping rates with Doppler applied
    # L5I
    f_l5i_d = L5_chipping_rate*(1. + f_d/sig_freq)
    f_l5i_dd = L5_chipping_rate*fd_rate/sig_freq
    # Neuman sequence
    f_nh_d = nh_chipping_rate*(1. + f_d/sig_freq)
    f_nh_dd = nh_chipping_rate*fd_rate/sig_freq
    # Databit
    # Data bit sequence
    f_db_d = L5_db_chipping_rate*(1. + f_d/sig_freq)
    f_db_dd = L5_db_chipping_rate*fd_rate/sig_freq
    # Calculate the L5I and nh code phase offsets
    l5i_init_code_phase = calcinitcodephase(L5_code_length,
                                            f_l5i_d, f_l5i_dd,
                                            f_s, code_start_idx)
    nh_init_code_phase = calcinitcodephase(nh10_code_length,
                                           f_nh_d, f_nh_dd,
                                           f_s, code_start_idx)
    db_init_code_phase = calcinitcodephase(db_code_length,
                                           f_db_d, f_db_dd,
                                           f_s, code_start_idx)
    # Allocate space for signal
    data = Array{Complex{Float64}}(undef, sample_num)
    isreplica = false
    noexp = false
    # Generate random databit vector
    databits = rand(0:1, db_code_length)
    return L5ISignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                     Tsys, CN0, ϕ, nADC, B, code_start_idx,
                     l5i_init_code_phase, nh_init_code_phase, db_init_code_phase,
                     t, data, include_carrier, include_adc, include_noise,
                     include_databits, f_l5i_d, f_l5i_dd,
                     f_nh_d, f_nh_dd, f_db_d, f_db_dd, sample_num,
                     isreplica, noexp, L5_chipping_rate,
                     sig_freq, L5_code_length, db_code_length, databits)
end


"""
    definesignal!(signal::L5ISignal;
                  prn=signal.prn, f_d=signal.f_d,
                  f_if=signal.f_if, fd_rate=signal.fd_rate,
                  Tsys=signal.Tsys, CN0=signal.CN0,
                  ϕ=signal.ϕ, nADC=signal.nADC,
                  B=signal.B,
                  include_carrier=signal.include_carrier,
                  include_adc=signal.include_adc,
                  include_noise=signal.include_noise,
                  code_start_idx=signal.code_start_idx,
                  isreplica=signal.isreplica,
                  noexp=signal.noexp,
                  include_databits=signal.include_databits,
                  sig_freq=signal.sig_freq)

Redefine properties of locally generated L5Q signal
based off its type, PRN, etc.

Sampling rate (f_s) and signal length (t_length)
are the only parameters that cannot be redefined.

"""
function definesignal!(signal::L5ISignal;
                       prn=signal.prn, f_d=signal.f_d,
                       f_if=signal.f_if, fd_rate=signal.fd_rate,
                       Tsys=signal.Tsys, CN0=signal.CN0,
                       ϕ=signal.ϕ, nADC=signal.nADC,
                       B=signal.B,
                       include_carrier=signal.include_carrier,
                       include_adc=signal.include_adc,
                       include_noise=signal.include_noise,
                       code_start_idx=signal.code_start_idx,
                       isreplica=signal.isreplica,
                       noexp=signal.noexp,
                       include_databits=signal.include_databits,
                       sig_freq=signal.sig_freq)
    # Calculate code chipping rates with Doppler applied
    # L5I
    f_l5i_d = L5_chipping_rate*(1. + f_d/sig_freq)
    f_l5i_dd = L5_chipping_rate*fd_rate/sig_freq
    # Neuman sequence
    f_nh_d = nh_chipping_rate*(1. + f_d/sig_freq)
    f_nh_dd = nh_chipping_rate*fd_rate/sig_freq
    # Databit
    # Data bit sequence
    f_db_d = L5_db_chipping_rate*(1. + f_d/sig_freq)
    f_db_dd = L5_db_chipping_rate*fd_rate/sig_freq
    # Calculate the L5I and nh code phase offsets
    l5i_init_code_phase = calcinitcodephase(L5_code_length,
                                            f_l5i_d, f_l5i_dd,
                                            signal.f_s, code_start_idx)
    nh_init_code_phase = calcinitcodephase(nh10_code_length,
                                           f_nh_d, f_nh_dd,
                                           signal.f_s, code_start_idx)
    db_init_code_phase = calcinitcodephase(signal.db_code_length,
                                           f_db_d, f_db_dd,
                                           signal.f_s, code_start_idx)
    # Store udated variables to "L5ISignal" struct
    signal.prn = prn
    signal.f_d = f_d
    signal.f_if = f_if
    signal.fd_rate = fd_rate
    signal.Tsys = Tsys
    signal.CN0 = CN0
    signal.ϕ = ϕ
    signal.nADC = nADC
    signal.B = B
    signal.include_carrier = include_carrier
    signal.include_adc = include_adc
    signal.include_noise = include_noise
    signal.f_l5i_d = f_l5i_d
    signal.f_l5i_dd = f_l5i_dd
    signal.f_nh_d = f_nh_d
    signal.f_nh_dd = f_nh_dd
    signal.f_db_d = f_db_d
    signal.f_db_dd = f_db_dd
    signal.l5i_init_code_phase = l5i_init_code_phase
    signal.nh_init_code_phase = nh_init_code_phase
    signal.db_init_code_phase = db_init_code_phase
    signal.code_start_idx = code_start_idx
    signal.isreplica = isreplica
    signal.noexp = noexp
    signal.include_databits = include_databits
    signal.sig_freq = sig_freq
    return signal
end
