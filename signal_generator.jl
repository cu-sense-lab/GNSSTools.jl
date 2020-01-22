using Random
using Base.Threads
include("constants.jl")
include("l5q_code_generator.jl")


"""
    L5QSignal()

Struct for holding L5Q GNSS signal properties for
signal generation.
"""
mutable struct L5QSignal{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,
                         A11,A12,A13,A14,A15,A16,A17,A18,
                         A19,A20,A21,A22,A23,A24,A25,A26,
                         A27,A28}
	type::A1
	prn::A2
	f_s::A3
	t_length::A4
	f_if::A5
	f_d::A6
	fd_rate::A7
	Tsys::A8
	CN0::A9
	ϕ::A10
	nADC::A11
	B::A12
	code_start_idx::A13
	l5q_init_code_phase::A14
	nh_init_code_phase::A15
	t::A16
	signal::A17
	include_carrier::A18
	include_adc::A19
	include_noise::A20
	include_neuman_code::A21
	f_l5q_d::A22
	f_l5q_dd::A23
	f_nh_d::A24
	f_nh_dd::A25
	sample_num::A26
	include_carrier_amplitude::A27
	isreplica::A28
end


"""
    calcinitcodephase(code_length, f_d, fd_rate, f_s,
                      code_start_idx)

Calculates the initial code phase of a given code
where f_d and fd_rate are the Doppler affected
code frequency and code frequency rate, respectively. 
"""
function calcinitcodephase(code_length, f_d, fd_rate, f_s,
                           code_start_idx)
	t₀ = (code_start_idx-1)/f_s
	init_phase = 1 - f_d*t₀ - 0.5*fd_rate*t₀^2
	if init_phase <= 0
		return code_length + 1 - abs(init_phase)%code_length
	else
		return init_phase%code_length
	end
end


"""
    definesignal(type::Val{:l5q}, prn, f_s, t_length;
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., ϕ=0., nADC=4, B=2.046e7,
                 include_carrier=true, include_adc=true,
                 include_noise=true, include_neuman_code=true,
                 code_start_idx=1, include_carrier_amplitude=true)

Define properties of locally generated L5Q signal
based off its type, PRN, etc.
"""
function definesignal(type::Val{:l5q}, prn, f_s, t_length;
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=4, B=2.046e7,
                      include_carrier=true, include_adc=true,
                      include_noise=true, include_neuman_code=true,
                      code_start_idx=1, include_carrier_amplitude=true)
	# Generate time vector
	t = Array(0:1/f_s:t_length-1/f_s)  # s
	sample_num = Int(f_s * t_length)
	# Calculate code chipping rates with Doppler applied
	# L5Q
	f_l5q_d = L5_chipping_rate*(1. + f_d/L5_freq)
	f_l5q_dd = L5_chipping_rate*fd_rate/L5_freq
	# Neuman sequence
	f_nh_d = nh_chipping_rate*(1. + f_d/L5_freq)
	f_nh_dd = nh_chipping_rate*fd_rate/L5_freq
	# Calculate the L5Q and nh code phase offsets
	l5q_init_code_phase = calcinitcodephase(L5_code_length,
                                            f_l5q_d, f_l5q_dd,
                                            f_s, code_start_idx)
	nh_init_code_phase = calcinitcodephase(nh_code_length,
                                           f_nh_d, f_nh_dd,
                                           f_s, code_start_idx)
	# Allocate space for signal
	signal = Array{Complex{Float64}}(undef, sample_num)
	isreplica = false
	return L5QSignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                     Tsys, CN0, ϕ, nADC, B, code_start_idx,
                     l5q_init_code_phase, nh_init_code_phase, t,
                     signal, include_carrier, include_adc,
                     include_noise, include_neuman_code,
                     f_l5q_d, f_l5q_dd,
                     f_nh_d, f_nh_dd, sample_num,
                     include_carrier_amplitude, isreplica)
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
                  include_neuman_code=signal.include_neuman_code,
                  code_start_idx=signal.code_start_idx,
                  include_carrier_amplitude=signal.include_carrier_amplitude)

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
                       include_neuman_code=signal.include_neuman_code,
                       code_start_idx=signal.code_start_idx,
                       include_carrier_amplitude=signal.include_carrier_amplitude,
                       isreplica=signal.isreplica)
	## Calculate code chipping rates with Doppler applied
	# L5Q
	f_l5q_d = L5_chipping_rate*(1. + f_d/L5_freq)
	f_l5q_dd = L5_chipping_rate*fd_rate/L5_freq
	# Neuman sequence
	f_nh_d = nh_chipping_rate*(1. + f_d/L5_freq)
	f_nh_dd = nh_chipping_rate*fd_rate/L5_freq
	# Calculate the L5Q and nh code phase offsets
	l5q_init_code_phase = calcinitcodephase(L5_code_length,
                                            f_l5q_d, f_l5q_dd,
                                            signal.f_s,
                                            code_start_idx)
	nh_init_code_phase = calcinitcodephase(nh_code_length,
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
	signal.include_neuman_code = signal.include_neuman_code
	signal.f_l5q_d = f_l5q_d
	signal.f_l5q_dd = f_l5q_dd
	signal.f_nh_d = f_nh_d
	signal.f_nh_dd = f_nh_dd
	signal.l5q_init_code_phase = l5q_init_code_phase
	signal.nh_init_code_phase = nh_init_code_phase
	signal.include_carrier_amplitude = include_carrier_amplitude
	signal.isreplica = isreplica
	return signal
end


"""
    calccodeidx(init_chip, f_code_d, f_code_dd, t, code_length)

Calculates the index in the codes for a given t.
"""
function calccodeidx(init_chip, f_code_d, f_code_dd,
                     t, code_length)
	chip = Int(floor(init_chip-1+f_code_d*t+0.5*f_code_dd*t^2)%code_length)+1
	if chip == 0
		return code_length
	else
		return chip
	end
end


"""
    generatesignal!(signal::L5QSignal,
                    isreplica::Val{false}=Val(signal.isreplica))

Generates local GNSS signal using paramters defined in a
L5QSignal struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

No need to specify `isreplica`. Set `isreplica` to `true` to
use alternate method, which ignores all `Bool` flags in `signal`.
"""
function generatesignal!(signal::L5QSignal,
	                     isreplica::Val{false}=Val(signal.isreplica))
	# Common parmeters used for entire signal
	prn = signal.prn
	Tsys = signal.Tsys
	CN0 = signal.CN0
	f_d = signal.f_d
	f_if = signal.f_if
	fd_rate = signal.fd_rate
	ϕ = signal.ϕ
	B = signal.B
	nADC = signal.nADC
	include_carrier = signal.include_carrier
	include_noise = signal.include_noise
	include_neuman_code = signal.include_neuman_code
	include_carrier_amplitude = signal.include_carrier_amplitude
	include_adc = signal.include_adc
	sigtype = eltype(signal.signal)
	adc_scale = 2^(nADC-1)-1
	@threads for i in 1:signal.sample_num
		@inbounds t = signal.t[i]
		# Get L5Q code value at t
		l5q = l5q_codes[prn][calccodeidx(signal.l5q_init_code_phase,
                                         signal.f_l5q_d, signal.f_l5q_dd,
                                         t, L5_code_length)]
		# Get Neuman code sequence value at t
		if include_neuman_code
			nh = nh20[calccodeidx(signal.nh_init_code_phase,
                                  signal.f_nh_d, signal.f_nh_dd,
                                  t, nh_code_length)]
		else
			nh = 0
		end
		if include_carrier & include_noise
			# Calculate code value with carrier and noise
			if include_carrier_amplitude
				@inbounds signal.signal[i] = ((xor(l5q, nh)*2-1)*sqrt(2*k*Tsys)*10^(CN0/20) *
			                                  exp((2π*(f_if + f_d + fd_rate*t)*t + ϕ)*1im) +
			                                  sqrt(k*B*Tsys)*randn(sigtype))
			else
				@inbounds signal.signal[i] = ((xor(l5q, nh)*2-1) *
			                                  exp((2π*(f_if + f_d + fd_rate*t)*t + ϕ)*1im) +
			                                  sqrt(k*B*Tsys)*randn(sigtype))
			end
		elseif include_carrier & ~include_noise
			# Calculate code value with carrier
			if include_carrier_amplitude
				@inbounds signal.signal[i] = ((xor(l5q, nh)*2-1)*sqrt(2*k*Tsys)*10^(CN0/20) *
			                                  exp((2π*(f_if + f_d + fd_rate*t)*t + ϕ)*1im))
			else
				@inbounds signal.signal[i] = ((xor(l5q, nh)*2-1) *
			                                  exp((2π*(f_if + f_d + fd_rate*t)*t + ϕ)*1im))
			end
		elseif ~include_carrier & include_noise
			# Calculate code value with carrier
			@inbounds signal.signal[i] = ((xor(l5q, nh)*2-1) +
			                              sqrt(k*B*Tsys)*randn(sigtype))
		else
			# Calculate code value only
			@inbounds signal.signal[i] = complex(float((xor(l5q, nh)*2-1)))
		end
	end
	# Quantize signal
	if include_adc
		signal.signal = round.(signal.signal.*adc_scale./sqrt(maximum(abs2.(signal.signal))))
	end
	return signal
end


"""
    generatesignal!(signal::L5QSignal,
                    isreplica::Val{true}=Val(signal.isreplica))

Generates local GNSS signal using paramters defined in a
L5QSignal struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.

This version is used only when `isreplica` is set to `true`
in `signal` and ignores all the `include_*` flags in `signal`,
except `include_neuman_code`. Exponential without the amplitude
is included automatically.
"""
function generatesignal!(signal::L5QSignal,
	                     isreplica::Val{true}=Val(signal.isreplica))
	# Common parmeters used for entire signal
	prn = signal.prn
	f_d = signal.f_d
	f_if = signal.f_if
	fd_rate = signal.fd_rate
	ϕ = signal.ϕ
	include_neuman_code = signal.include_neuman_code
	@threads for i in 1:signal.sample_num
		@inbounds t = signal.t[i]
		# Get L5Q code value at t
		l5q = l5q_codes[prn][calccodeidx(signal.l5q_init_code_phase,
                                         signal.f_l5q_d, signal.f_l5q_dd,
                                         t, L5_code_length)]
		# Get Neuman code sequence value at t
		if include_neuman_code
			nh = nh20[calccodeidx(signal.nh_init_code_phase,
                                  signal.f_nh_d, signal.f_nh_dd,
                                  t, nh_code_length)]
		else
			nh = 0
		end
		@inbounds signal.signal[i] = ((xor(l5q, nh)*2-1) *
			                          exp((2π*(f_if + f_d + fd_rate*t)*t + ϕ)*1im))
	end
	return signal
end
