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
    if ~signal.signal_type.include_I
        B = signal.signal_type.B_Q
    elseif ~signal.signal_type.include_Q
        B = signal.signal_type.B_I
    else
        B = max(signal.signal_type.B_I, signal.signal_type.B_Q)
    end
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
    upsample_factor = denominator(signal.code_start_idx)
    Δt = signal.t[2] - signal.t[1]
    dΔt = Δt/upsample_factor
    @threads for i in 1:N
        @inbounds t = signal.t[i]
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
            ϕ_noise = phase_noise[i]
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ+ϕ_noise) +
            #                            noise_amp * thermal_noise[i]
            @inbounds signal.data[i] = carrier_amp * cis_sum * cis(ϕ_noise) +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and noise
            ϕ = get_ϕ(t)
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ) +
            #                            noise_amp * thermal_noise[i]
            @inbounds signal.data[i] = carrier_amp * cis_sum +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & ~include_thermal_noise & include_phase_noise
            # Calculate code value with noise and no carrier
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ+ϕ_noise)
            @inbounds signal.data[i] = carrier_amp * cis_sum * cis(ϕ_noise)
        elseif include_carrier & ~include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and no noise
            ϕ = get_ϕ(t)
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ)
            @inbounds signal.data[i] = carrier_amp * cis_sum
        else
            # Calculate code value only
            # @inbounds signal.data[i] = complex(float(code_val))
            @inbounds signal.data[i] = cis_sum
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
    if ~signal.signal_type.include_I
        B = signal.signal_type.B_Q
    elseif ~signal.signal_type.include_Q
        B = signal.signal_type.B_I
    else
        B = max(signal.signal_type.B_I, signal.signal_type.B_Q)
    end
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
    upsample_factor = denominator(signal.code_start_idx)
    Δt = signal.t[2] - signal.t[1]
    dΔt = Δt/upsample_factor
    @threads for i in 1:N
        @inbounds t = signal.t[i]
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
            ϕ_noise = phase_noise[i]
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ+ϕ_noise) +
            #                            noise_amp * thermal_noise[i]
            @inbounds signal.data[i] = carrier_amp * cis_sum * cis(ϕ_noise) +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and noise
            ϕ = get_ϕ(t)
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ) +
            #                            noise_amp * thermal_noise[i]
            @inbounds signal.data[i] = carrier_amp * cis_sum +
                                       noise_amp * thermal_noise[i]
        elseif include_carrier & ~include_thermal_noise & include_phase_noise
            # Calculate code value with noise and no carrier
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ+ϕ_noise)
            @inbounds signal.data[i] = carrier_amp * cis_sum * cis(ϕ_noise)
        elseif include_carrier & ~include_thermal_noise & ~include_phase_noise
            # Calculate code value with carrier and no noise
            ϕ = get_ϕ(t)
            # @inbounds signal.data[i] = code_val * carrier_amp * cis(ϕ)
            @inbounds signal.data[i] = carrier_amp * cis_sum
        else
            # Calculate code value only
            # @inbounds signal.data[i] = complex(float(code_val))
            @inbounds signal.data[i] = cis_sum
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
    generatesignal!(signal::ReplicaSignal, isreplica::Bool)


Generates local GNSS signal using paramters defined in a
`ReplicaSignal` struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

This version is used only when `isreplica` is set to `true`
in `signal` and ignores all the `include_*` flags in `signal`.
Exponential without the amplitude is included automatically.


Required Arguments:

- `signal::ReplicaSignal`: struct containing signal data and parameters
- `isreplica::Bool`: typically set to `true` and clarifies that `signal` is
                     being used for processing only with no noise added


Modifies in Place and Returns:

- `signal::ReplicaSignal`: generated signal data stored in `signal.data`
"""
function generatesignal!(signal::ReplicaSignal, isreplica::Bool)
    # Common parmeters used for entire signal
    prn = signal.prn
    f_d = signal.f_d
    f_if = signal.f_if
    fd_rate = signal.fd_rate
    ϕ = signal.phi
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
