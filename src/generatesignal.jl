"""
    generatesignal!(signal::ReplicaSignal, t_length=signal.t_length;
                    doppler_curve=missing, doppler_t=missing)


Generates local GNSS signal using paramters defined in a `ReplicaSignal`
struct.

Generates a signal with carrier, ADC quantization, noise, and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to use alternate
method, which ignores all `Bool` flags in `signal`.

Specify `t_length` to be ≲ to `replica.t_length` for different signal generation
lengths. **NOTE:** The size of the `replica.data` array will still be the same
size. You will need to keep track of the `t_length` you passed to
`generatesignal!`.


Required Arguments:

- `signal::ReplicaSignal`: struct containing signal data and parameters


Optional Arguments:

- `t_length::Float64`: signal generation length in seconds
                       `(default = signal.t_length)`
- `doppler_curve::Vector`: Doppler frequency curve in Hz to use instead of
                           `signal.f_d` and `signal.fd_rate` parameters
                           `(default = missing)`
- `doppler_t::Vector`: time vector in seconds corresponding to `doppler_curve`
                       `(default = missing)`


Modifies in Place:

- `signal::ReplicaSignal`: generated signal data stored in `signal.data`
"""
function generatesignal!(signal::ReplicaSignal, t_length=signal.t_length;
                         doppler_curve=missing, doppler_t=missing)
    # Common parmeters used for entire signal
    prn = signal.prn
    Tsys = signal.Tsys
    CN0 = signal.CN0
    f_s = signal.f_s
    f_if = signal.f_if
    f_d = signal.f_d
    fd_rate = signal.fd_rate
    ϕ_init = signal.phi
    if ismissing(doppler_curve)
        get_ϕ(t) = 2π*(f_if + f_d + 0.5*fd_rate*t)*t + ϕ_init
        get_code_val(t) = calc_code_val(signal, t)
    else
        get_code_val, get_ϕ = get_chips_and_ϕ(signal, doppler_curve, doppler_t)
    end
    generatesignal!(signal, t_length, get_code_val, get_ϕ)
end


"""
    generatesignal!(signal::ReplicaSignal, t_length, get_code_val, get_ϕ)


Generates local GNSS signal using paramters defined in a `ReplicaSignal`
struct.

Generates a signal with carrier, ADC quantization, noise, and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to use alternate
method, which ignores all `Bool` flags in `signal`.

Specify `t_length` to be ≲ to `replica.t_length` for different signal generation
lengths. **NOTE:** The size of the `replica.data` array will still be the same
size. You will need to keep track of the `t_length` you passed to
`generatesignal!`.

`get_code_val` and `get_ϕ` are both only functions of `t`, where `t` is in
seconds. They return the current code value for all codes and carrier phase,
respectively.


Required Arguments:

- `signal::ReplicaSignal`: struct containing signal data and parameters
- `t_length::Float64`: signal generation length in seconds
- `get_code_val`: a function that calculates the current code value for a given
                  `t`
    * Required Arguments:
        - `t::Float64`: current time in seconds
- `get_ϕ`: a function that calculates the current carrier phase for a given `t`
    * Required Arguments:
        - `t::Float64`: current time in seconds


Modifies in Place and Returns:

- `signal::ReplicaSignal`: generated signal data stored in `signal.data`
"""
function generatesignal!(signal::ReplicaSignal, t_length, get_code_val, get_ϕ)
    # Common parmeters used for entire signal
    prn = signal.prn
    Tsys = signal.Tsys
    CN0 = signal.CN0
    f_s = signal.f_s
    start_t = signal.start_t
    if ~signal.signal_type.include_I && signal.signal_type.include_Q
        B = signal.signal_type.B_Q
    elseif signal.signal_type.include_I && ~signal.signal_type.include_Q
        B = signal.signal_type.B_I
    else  # assume both channels are included and determine highest bandwidth
        B = max(signal.signal_type.B_I, signal.signal_type.B_Q)
    end
    nADC = signal.nADC
    include_carrier = signal.include_carrier
    include_thermal_noise = signal.include_thermal_noise
    include_phase_noise = signal.include_phase_noise
    include_adc = signal.include_adc
    thermal_noise = signal.thermal_noise
    phase_noise = signal.phase_noise
    adc_scale = 2^(nADC-1)
    carrier_amp = sqrt(2*k*Tsys)*10^(CN0/20)  # for sinusoid
    noise_amp = sqrt(k*B*Tsys)
    N = floor(Int, t_length*f_s)
    upsample_factor = denominator(signal.code_start_idx)
    Δt = 1/f_s
    dΔt = Δt/upsample_factor
    @threads for i in 1:N
        t = calc_t_at_i(i, start_t, f_s)
        # Generate code value for given signal type
        cis_sum = complex(0.)
        if include_carrier
            for j in 0:(upsample_factor-1)
                code_val, code_ϕ = get_code_val(t+j*dΔt)
                ϕ = get_ϕ(t+j*dΔt)
                cis_sum += cis(ϕ+code_ϕ)
                # cis_sum += real(code_val)*cos(ϕ) + imag(code_val)*sin(ϕ)*1im
            end
        else
            for j in 0:(upsample_factor-1)
                code_val, code_ϕ = get_code_val(t+j*dΔt)
                cis_sum += cis(code_ϕ)
                # cis_sum += real(code_val) + imag(code_val)*1im
            end
        end
        cis_sum = cis_sum/upsample_factor
        if include_carrier & include_thermal_noise & include_phase_noise
            # Calculate code value with carrier, thermal and phase noise
            ϕ_noise = real(phase_noise[i])
            @inbounds signal.data[i] = carrier_amp * cis_sum * cis(ϕ_noise) +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and noise
            @inbounds signal.data[i] = carrier_amp * cis_sum +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & ~include_thermal_noise & include_phase_noise
            # Calculate code value with noise and no carrier
            ϕ_noise = real(phase_noise[i])
            @inbounds signal.data[i] = carrier_amp * cis_sum * cis(ϕ_noise)
        elseif include_carrier & ~include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and no noise
            @inbounds signal.data[i] = carrier_amp * cis_sum
        else
            # Calculate code value only
            @inbounds signal.data[i] = cis_sum
        end
    end
    # Quantize signal
    if include_adc
        sigmax = sqrt(maximum(abs2, signal.data))
        @threads for i in 1:signal.sample_num
            @inbounds signal.data[i] = round(signal.data[i]*adc_scale/sigmax)
            # Check that the max positive value for I/Q channels is 
            # (2^(nADC - 1) - 1)
            if real(signal.data[i]) == adc_scale
                signal.data[i] = signal.data[i] - 1 
            end
            if imag(signal.data[i]) == adc_scale
                signal.data[i] = signal.data[i] - 1im
            end
        end
    end
    return signal
end


"""
    generatesignal!(signal::ReplicaSignals, t_length=signal.t_length;
                    doppler_curve=missing, doppler_t=missing)


Generates local GNSS signal using paramters defined in a `ReplicaSignals`
struct.

Generates a signal with carrier, ADC quantization, noise, and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to use alternate
method, which ignores all `Bool` flags in `signal`.

Specify `t_length` to be ≲ to `replica.t_length` for different signal generation
lengths. **NOTE:** The size of the `replica.data` array will still be the same
size. You will need to keep track of the `t_length` you passed to
`generatesignal!`.


Required Arguments:

- `signal::ReplicaSignals`: struct containing signal data and parameters


Optional Arguments:

- `t_length::Float64`: signal generation length in seconds
                       `(default = signal.t_length)`
- `doppler_curve::Vector`: Doppler frequency curve in Hz to use instead of
                           `signal.f_d` and `signal.fd_rate` parameters
                           `(default = missing)`
- `doppler_t::Vector`: time vector in seconds corresponding to `doppler_curve`
                       `(default = missing)`


Modifies in Place:

- `signal::ReplicaSignal`: generated signal data stored in `signal.data`
"""
function generatesignal!(signal::ReplicaSignals, t_length=signal.t_length;
                         doppler_curve=missing, doppler_t=missing)
    # Common parmeters used for entire signal
    prn = signal.prn
    Tsys = signal.Tsys
    CN0 = signal.CN0
    f_s = signal.f_s
    f_if = signal.f_if
    f_d = signal.f_d
    fd_rate = signal.fd_rate
    ϕ_init = signal.phi
    if ismissing(doppler_curve)
        get_ϕ(t) = 2π*(f_if + f_d + 0.5*fd_rate*t)*t + ϕ_init
        get_code_val(t) = calc_code_val(signal, t)
    else
        get_code_val, get_ϕ = get_chips_and_ϕ(signal, doppler_curve, doppler_t)
    end
    generatesignal!(signal, t_length, get_code_val, get_ϕ)
end


"""
    generatesignal!(signal::ReplicaSignals, t_length, get_code_val, get_ϕ)


Generates local GNSS signal using paramters defined in a `ReplicaSignals`
struct.

Generates a signal with carrier, ADC quantization, noise, and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to use alternate
method, which ignores all `Bool` flags in `signal`.

Specify `t_length` to be ≲ to `replica.t_length` for different signal generation
lengths. **NOTE:** The size of the `replica.data` array will still be the same
size. You will need to keep track of the `t_length` you passed to
`generatesignal!`.

`get_code_val` and `get_ϕ` are both only functions of `t`, where `t` is in
seconds. They return the current code value for all codes and carrier phase,
respectively.


Required Arguments:

- `signal::ReplicaSignals`: struct containing signal data and parameters
- `t_length::Float64`: signal generation length in seconds
- `get_code_val`: a function that calculates the current code value for a given
                  `t`
    * Required Arguments:
        - `t::Float64`: current time in seconds
- `get_ϕ`: a function that calculates the current carrier phase for a given `t`
    * Required Arguments:
        - `t::Float64`: current time in seconds


Modifies in Place and Returns:

- `signal::ReplicaSignal`: generated signal data stored in `signal.data`
"""
function generatesignal!(signal::ReplicaSignals, t_length, get_code_val, get_ϕ)
    # Common parmeters used for entire signal
    prn = signal.prn
    Tsys = signal.Tsys
    CN0 = signal.CN0
    f_s = signal.f_s
    start_t = signal.start_t
    B = signal.B
    nADC = signal.nADC
    include_carrier = signal.include_carrier
    include_thermal_noise = signal.include_thermal_noise
    include_phase_noise = signal.include_phase_noise
    include_adc = signal.include_adc
    thermal_noise = signal.thermal_noise
    phase_noise = signal.phase_noise
    adc_scale = 2^(nADC-1)
    carrier_amp = sqrt(2*k*Tsys)*10^(CN0/20)
    noise_amp = sqrt(k*B*Tsys)
    N = floor(Int, t_length*f_s)
    upsample_factor = denominator(signal.code_start_idx)
    Δt = 1/f_s
    dΔt = Δt/upsample_factor
    @threads for i in 1:N
        t = calc_t_at_i(i, start_t, f_s)
        # Generate code value for given signal type
        cis_sum = complex(0.)
        if include_carrier
            for j in 0:(upsample_factor-1)
                code_val, code_ϕ = get_code_val(t+dΔt)
                ϕ = get_ϕ(t+dΔt)
                cis_sum += cis(ϕ+code_ϕ)
            end
        else
            for j in 0:(upsample_factor-1)
                code_val, code_ϕ = get_code_val(t+dΔt)
                cis_sum += cis(code_ϕ)
            end
        end
        cis_sum = cis_sum/upsample_factor
        if include_carrier & include_thermal_noise & include_phase_noise
            # Calculate code value with carrier, thermal and phase noise
            ϕ = get_ϕ(t)
            ϕ_noise = real(phase_noise[i])
            @inbounds signal.data[i] = carrier_amp * cis_sum * cis(ϕ_noise) +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and noise
            # ϕ = get_ϕ(t)
            @inbounds signal.data[i] = carrier_amp * cis_sum +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & ~include_thermal_noise & include_phase_noise
            # Calculate code value with noise and no carrier
            ϕ_noise = real(phase_noise[i])
            @inbounds signal.data[i] = carrier_amp * cis_sum * cis(ϕ_noise)
        elseif include_carrier & ~include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and no noise
            ϕ = get_ϕ(t)
            @inbounds signal.data[i] = carrier_amp * cis_sum
        else
            # Calculate code value only
            @inbounds signal.data[i] = cis_sum
        end
    end
    # Quantize signal
    if include_adc
        sigmax = sqrt(maximum(abs2, signal.data))
        # for i in 1:signal.sample_num
        @threads for i in 1:signal.sample_num
            @inbounds signal.data[i] = round(signal.data[i]*adc_scale/sigmax)
        end
    end
    return signal
end


"""
generatereplica!(signal::ReplicaSignal)


Generates local GNSS signal using paramters defined in a
`ReplicaSignal` struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

This version is used only when for when this signal struct is
used as a known replica for use in course/fine acquisition and 
tracking. When this function is used instead of `generatesignal!`,
all the `include_*` flags in `signal` are ignored. The bool flag,
`noexp` in `signal`, Exponential without the amplitude is included 
automatically.


Required Arguments:

- `signal::ReplicaSignal`: struct containing signal data and parameters


Modifies in Place and Returns:

- `signal::ReplicaSignal`: generated signal data stored in `signal.data`
"""
function generatereplica!(signal::ReplicaSignal)
    # Common parmeters used for entire signal
    prn = signal.prn
    f_d = signal.f_d
    f_if = signal.f_if
    f_s = signal.f_s
    start_t = signal.start_t
    fd_rate = signal.fd_rate
    ϕ = signal.phi
    noexp = signal.noexp
    signal.isreplica = true
    @threads for i in 1:signal.sample_num
        t = calc_t_at_i(i, start_t, f_s)
        # Generate code value for given signal type
        code_val, code_ϕ = calc_code_val(signal, t)
        if signal.include_carrier
            # Doppler effected code is modulated onto carrier.
            # Carrier has amplitude of 1. It is not rescaled.
            ϕₜ = 2π*(f_if + f_d + 0.5*fd_rate*t)*t + ϕ
            @inbounds signal.data[i] = cis(ϕₜ + code_ϕ)
            # @inbounds signal.data[i] = real(code_val)*cos(ϕₜ) + imag(code_val)*sin(ϕₜ)*1im
        else
            # Produce the code only with no carrier modulation
            # Code is still generated with Doppler effects
            @inbounds signal.data[i] = cis(code_ϕ)
            # @inbounds signal.data[i] = real(code_val) + imag(code_val)*1im
        end
    end
    signal.isreplica = false
    return signal
end
