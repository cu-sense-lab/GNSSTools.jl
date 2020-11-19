"""
    generatesignal!(signal::ReplicaSignals,
                    t_length=signal.t_length;
                    doppler_curve=missing, doppler_t=missing,
                    message="Generating signal...")

Generates local GNSS signal using paramters defined in a
`ReplicaSignals` struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to
use alternate method, which ignores all `Bool` flags in `signal`.

Specify `t_length` to be ≲ to `replica.t_length` for different
signal generation lengths. **NOTE:** The size of the `replica.data`
array will still be the same size. You will need to keep track
of the `t_length` you passed to `generatesignal!`.
"""
function generatesignal!(signal::ReplicaSignals,
                         t_length=signal.t_length;
                         doppler_curve=missing, doppler_t=missing,
                         message="Generating signal...")
    # Common parmeters used for entire signal
    prn = signal.prn
    Tsys = signal.Tsys
    CN0 = signal.CN0
    f_s = signal.f_s
    f_if = signal.f_if
    f_d = signal.f_d
    fd_rate = signal.fd_rate
    ϕ_init = signal.ϕ
    if ismissing(doppler_curve)
        get_ϕ(t) = 2π*(f_if + f_d + 0.5*fd_rate*t)*t + ϕ_init
        get_code_val(t) = calc_code_val(signal, t)
    else
        get_code_val, get_ϕ = get_chips_and_ϕ(signal, doppler_curve;
                                              doppler_t=doppler_t)
    end
    generatesignal!(signal, t_length, get_code_val, get_ϕ)
end


"""
    generatesignal!(signal::ReplicaSignals, t_length, get_code_val, get_ϕ)

Generates local GNSS signal using paramters defined in a
`ReplicaSignals` struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to
use alternate method, which ignores all `Bool` flags in `signal`.

Specify `t_length` to be ≲ to `replica.t_length` for different
signal generation lengths. **NOTE:** The size of the `replica.data`
array will still be the same size. You will need to keep track
of the `t_length` you passed to `generatesignal!`.
"""
function generatesignal!(signal::ReplicaSignals, t_length, get_code_val, get_ϕ)
    # Common parmeters used for entire signal
    prn = signal.prn
    Tsys = signal.Tsys
    CN0 = signal.CN0
    f_s = signal.f_s
    B = signal.signal_type.B
    nADC = signal.nADC
    include_carrier = signal.include_carrier
    include_thermal_noise = signal.include_thermal_noise
    include_phase_noise = signal.include_phase_noise
    include_adc = signal.include_adc
    thermal_noise = signal.thermal_noise
    phase_noise = signal.phase_noise
    sigtype = eltype(signal.data)
    adc_scale = 2^(nADC-1)-1
    carrier_amp = sqrt(2*k*Tsys)*10^(CN0/20)
    noise_amp = sqrt(k*B*Tsys)
    N = floor(Int, t_length*f_s)
    @threads for i in 1:N
        @inbounds t = signal.t[i]
        # Generate code value for given signal type
        code_val, code_ϕ = get_code_val(t)
        if include_carrier & include_thermal_noise & include_phase_noise
            # Calculate code value with carrier, thermal and phase noise
            ϕ = get_ϕ(t)
            ϕ_noise = phase_noise[i]
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ+ϕ_noise) +
            #                            noise_amp * thermal_noise[i]
            @inbounds signal.data[i] = carrier_amp * cis(ϕ+ϕ_noise+code_ϕ) +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and noise
            ϕ = get_ϕ(t)
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ) +
            #                            noise_amp * thermal_noise[i]
            @inbounds signal.data[i] = carrier_amp * cis(ϕ+code_ϕ) +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & ~include_thermal_noise & include_phase_noise
            # Calculate code value with noise and no carrier
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ+ϕ_noise)
            @inbounds signal.data[i] = carrier_amp * cis(ϕ+ϕ_noise+code_ϕ)
        elseif include_carrier & ~include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and no noise
            ϕ = get_ϕ(t)
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ)
            @inbounds signal.data[i] = carrier_amp * cis(ϕ+code_ϕ)
        else
            # Calculate code value only
            # @inbounds signal.data[i] = complex(float(code_val))
            @inbounds signal.data[i] = cis(code_ϕ)
        end
    end
    # Quantize signal
    if include_adc
        sigmax = sqrt(maximum(abs2, signal.data))
        @threads for i in 1:signal.sample_num
            @inbounds signal.data[i] = round(signal.data[i]*adc_scale/sigmax)
        end
    end
    return signal
end


"""
    generatesignal!(signal::ReplicaSignals, isreplica::Bool)

Generates local GNSS signal using paramters defined in a
`ReplicaSignals` struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

This version is used only when `isreplica` is set to `true`
in `signal` and ignores all the `include_*` flags in `signal`.
Exponential without the amplitude is included automatically.
"""
function generatesignal!(signal::ReplicaSignals, isreplica::Bool)
    # Common parmeters used for entire signal
    prn = signal.prn
    f_d = signal.f_d
    f_if = signal.f_if
    fd_rate = signal.fd_rate
    ϕ = signal.ϕ
    noexp = signal.noexp
    @threads for i in 1:signal.sample_num
        @inbounds t = signal.t[i]
        # Generate code value for given signal type
        code_val, code_ϕ = calc_code_val(signal, t)
        if noexp
            # @inbounds signal.data[i] = complex(float(code_val))
            @inbounds signal.data[i] = cis(code_ϕ)
        else
            # @inbounds signal.data[i] = (code_val *
            #                             cis(2π*(f_if + f_d + 0.5*fd_rate*t)*t + ϕ))
            @inbounds signal.data[i] = cis(2π*(f_if + f_d + 0.5*fd_rate*t)*t + ϕ + code_ϕ)
        end
    end
    return signal
end


#------------------------------------------------------------------------------
#                             OLD IMPLEMENTATON
#------------------------------------------------------------------------------


"""
    generatesignal!(signal::ReplicaSignal,
                    t_length=signal.t_length;
                    doppler_curve=missing, doppler_t=missing,
                    message="Generating signal...")

Generates local GNSS signal using paramters defined in a
`ReplicaSignal` struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to
use alternate method, which ignores all `Bool` flags in `signal`.

Specify `t_length` to be ≲ to `replica.t_length` for different
signal generation lengths. **NOTE:** The size of the `replica.data`
array will still be the same size. You will need to keep track
of the `t_length` you passed to `generatesignal!`.
"""
function generatesignal!(signal::ReplicaSignal,
                         t_length=signal.t_length;
                         doppler_curve=missing, doppler_t=missing,
                         message="Generating signal...")
    # Common parmeters used for entire signal
    prn = signal.prn
    Tsys = signal.Tsys
    CN0 = signal.CN0
    f_s = signal.f_s
    f_if = signal.f_if
    f_d = signal.f_d
    fd_rate = signal.fd_rate
    ϕ_init = signal.ϕ
    if ismissing(doppler_curve)
        get_ϕ(t) = 2π*(f_if + f_d + 0.5*fd_rate*t)*t + ϕ_init
        get_code_val(t) = calc_code_val(signal, t)
    else
        get_code_val, get_ϕ = get_chips_and_ϕ(signal, doppler_curve;
                                              doppler_t=doppler_t)
    end
    N = floor(Int, t_length*f_s)
    sigtype = eltype(signal.data)
    thermal_noise = randn(sigtype, N)
    generatesignal!(signal, t_length, get_code_val, get_ϕ, thermal_noise)
end


"""
    generatesignal!(signal::ReplicaSignal, t_length, get_code_val, get_ϕ, thermal_noise)

Generates local GNSS signal using paramters defined in a
`ReplicaSignal` struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to
use alternate method, which ignores all `Bool` flags in `signal`.

Specify `t_length` to be ≲ to `replica.t_length` for different
signal generation lengths. **NOTE:** The size of the `replica.data`
array will still be the same size. You will need to keep track
of the `t_length` you passed to `generatesignal!`.
"""
function generatesignal!(signal::ReplicaSignal, t_length, get_code_val, get_ϕ, thermal_noise)
    # Common parmeters used for entire signal
    prn = signal.prn
    Tsys = signal.Tsys
    CN0 = signal.CN0
    f_s = signal.f_s
    B = signal.B
    nADC = signal.nADC
    include_carrier = signal.include_carrier
    include_noise = signal.include_noise
    include_adc = signal.include_adc
    sigtype = eltype(signal.data)
    adc_scale = 2^(nADC-1)-1
    carrier_amp = sqrt(2*k*Tsys)*10^(CN0/20)
    noise_amp = sqrt(k*B*Tsys)
    N = floor(Int, t_length*f_s)
    # p = Progress(N, 15, message)
    @threads for i in 1:N
        @inbounds t = signal.t[i]
        # Generate code value for given signal type
        code_val = get_code_val(t)
        if include_carrier & include_noise
            # Calculate code value with carrier and noise
            ϕ = get_ϕ(t)
            @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ) +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & ~include_noise
            # Calculate code value with carrier and no noise
            ϕ = get_ϕ(t)
            @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ)
        elseif ~include_carrier & include_noise
            # Calculate code value with noise and no carrier
            @inbounds signal.data[i] = code_val + noise_amp * thermal_noise[i]
        else
            # Calculate code value only
            @inbounds signal.data[i] = complex(float(code_val))
        end
        # next!(p)
    end
    # Quantize signal
    if include_adc
        sigmax = sqrt(maximum(abs2, signal.data))
        @threads for i in 1:signal.sample_num
            @inbounds signal.data[i] = round(signal.data[i]*adc_scale/sigmax)
        end
    end
    return signal
end


"""
    generatesignal!(signal::ReplicaSignal, isreplica::Bool)

Generates local GNSS signal using paramters defined in a
`ReplicaSignal` struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

This version is used only when `isreplica` is set to `true`
in `signal` and ignores all the `include_*` flags in `signal`.
Exponential without the amplitude is included automatically.
"""
function generatesignal!(signal::ReplicaSignal, isreplica::Bool)
    # Common parmeters used for entire signal
    prn = signal.prn
    f_d = signal.f_d
    f_if = signal.f_if
    fd_rate = signal.fd_rate
    ϕ = signal.ϕ
    noexp = signal.noexp
    @threads for i in 1:signal.sample_num
        @inbounds t = signal.t[i]
        # Generate code value for given signal type
        code_val = calc_code_val(signal, t)
        if noexp
            @inbounds signal.data[i] = complex(float(code_val))
        else
            @inbounds signal.data[i] = (code_val *
                                        exp((2π*(f_if + f_d + 0.5*fd_rate*t)*t + ϕ)*1im))
        end
    end
    return signal
end
