include("constants.jl")
include("l5q_code_generator.jl")


"""
L5QSignal()

Struct for holding L5Q GNSS signal properties for
signal generation.
"""
mutable struct L5QSignal{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,
                         A11,A12,A13,A14,A15,A16,A17,A18,
                         A19,A20,A21,A22,A23,A24,A25,A26}
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
definesignal(type, prn, f_s, t_length;
             f_if, f_d, fd_rate, Tsys,
             CN0, ϕ, nADC, B, include_carrier,
             include_adc, include_noise,
             include_neuman_code, code_start_idx,
             fds)

Define properties of locally generated L5Q signal
based off its type, PRN, etc.
"""
function definesignal(type::Val{:l5q}, prn, f_s, t_length;
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=4, B=2.046e7,
                      include_carrier=true, include_adc=true,
                      include_noise=true, include_neuman_code=true,
                      code_start_idx=1)
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
	return L5QSignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                     Tsys, CN0, ϕ, nADC, B, code_start_idx,
                     l5q_init_code_phase, nh_init_code_phase, t,
                     signal, include_carrier, include_adc,
                     include_noise, include_neuman_code,
                     f_l5q_d, f_l5q_dd,
                     f_nh_d, f_nh_dd, sample_num)
end


"""
definesignal!(l5qsignal; prn, f_if, f_d, fd_rate,
              Tsys, CN0, ϕ, nADC, B, include_carrier,
              include_adc, include_noise,
              include_neuman_code, code_start_idx, fds)

Redefine properties of locally generated L5Q signal
based off its type, PRN, etc. 

Sampling rate (f_s) and signal length (t_length)
are the only parameters that cannot be redefined.

"""
function definesignal!(l5qsignal::L5QSignal;
                       prn=l5qsignal.prn, f_d=l5qsignal.f_d,
                       f_if=l5qsignal.f_if, fd_rate=l5qsignal.fd_rate,
                       Tsys=l5qsignal.Tsys, CN0=l5qsignal.CN0,
                       ϕ=l5qsignal.ϕ, nADC=l5qsignal.nADC,
                       B=l5qsignal.B,
                       include_carrier=l5qsignal.include_carrier,
                       include_adc=l5qsignal.include_adc,
                       include_noise=l5qsignal.include_noise,
                       include_neuman_code=l5qsignal.include_neuman_code,
                       code_start_idx=l5qsignal.code_start_idx)
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
                                            l5qsignal.f_s,
                                            code_start_idx)
	nh_init_code_phase = calcinitcodephase(nh_code_length,
                                           f_nh_d, f_nh_dd,
                                           l5qsignal.f_s,
                                           code_start_idx)
	# Store udated variables to "l5qsignal" struct
	l5qsignal.prn = prn
	l5qsignal.f_d = f_d
	l5qsignal.f_if = f_if
	l5qsignal.fd_rate = fd_rate
	l5qsignal.Tsys = Tsys
	l5qsignal.CN0 = CN0
	l5qsignal.ϕ = ϕ
	l5qsignal.nADC = nADC
	l5qsignal.B = B
	l5qsignal.include_carrier = include_carrier
	l5qsignal.include_adc = include_adc
	l5qsignal.include_noise = include_noise
	l5qsignal.include_neuman_code = l5qsignal.include_neuman_code
	l5qsignal.f_l5q_d = f_l5q_d
	l5qsignal.f_l5q_dd = f_l5q_dd
	l5qsignal.f_nh_d = f_nh_d
	l5qsignal.f_nh_dd = f_nh_dd
	l5qsignal.l5q_init_code_phase = l5q_init_code_phase
	l5qsignal.nh_init_code_phase = nh_init_code_phase
	return l5qsignal
end


"""
calccodeidx(init_chip, f_code_d, f_code_d, t, code_length)

Calculates the index in the codes for a given t.
"""
function calccodeidx(init_chip, f_code_d, f_code_dd,
                     t, code_length)
	return Int(floor(init_chip+f_code_d*t+0.5*f_code_dd*t^2)%code_length)+1
end


"""
generatesignal!(l5qsignal)

Generates local GNSS signal using paramters defined in a
L5QSignal struct.

Generates a signal with carrier, ADC quantization, noise,
and Neuman sequence.
"""
function generatesignal!(l5qsignal::L5QSignal)#,
                         # include_carrier::Val{true}=Val(l5qsignal.include_carrier),
                         # include_noise::Val{true}=Val(l5qsignal.include_noise),
                         # include_neuman_code::Val{true}=Val(l5qsignal.include_neuman_code))
	# Parmeters used for entire signal
	prn = l5qsignal.prn
	Tsys = l5qsignal.Tsys
	CN0 = l5qsignal.CN0
	f_d = l5qsignal.f_d
	f_if = l5qsignal.f_if
	fd_rate = l5qsignal.fd_rate
	ϕ = l5qsignal.ϕ
	B = l5qsignal.B
	nADC = l5qsignal.nADC
	include_carrier = l5qsignal.include_carrier
	include_noise = l5qsignal.include_noise
	include_neuman_code = l5qsignal.include_neuman_code
	sigtype = eltype(l5qsignal.signal)
	for i in 1:l5qsignal.sample_num
		t = l5qsignal.t[i]
		# Get L5Q code value at t
		l5q = l5q_codes[prn][calccodeidx(l5qsignal.l5q_init_code_phase,
                                         l5qsignal.f_l5q_d, l5qsignal.f_l5q_dd,
                                         t, L5_code_length)]
		# Get Neuman code sequence value at t
		if include_neuman_code
			nh = nh20[calccodeidx(l5qsignal.nh_init_code_phase,
                                  l5qsignal.f_nh_d, l5qsignal.f_nh_dd,
                                  t, nh_code_length)]
		else
			nh = 0
		end
		if include_carrier & include_noise
			# Calculate code value with carrier and noise
			sig = ((xor(l5q, nh)*2-1)*sqrt(2*k*Tsys)*10^(CN0/20) *
			       exp((2π*(f_if + f_d + fd_rate*t)*t + ϕ)*1im) +
			       sqrt(k*B*Tsys)*randn(sigtype))
		elseif include_carrier & ~include_noise
			# Calculate code value with carrier
			sig = ((xor(l5q, nh)*2-1)*sqrt(2*k*Tsys)*10^(CN0/20) *
			       exp((2π*(f_if + f_d + fd_rate*t)*t + ϕ)*1im))
		elseif ~include_carrier & include_noise
			# Calculate code value with carrier
			sig = ((xor(l5q, nh)*2-1) +
			       sqrt(k*B*Tsys)*randn(sigtype))
		else
			# Calculate code value only
			sig = complex(float((xor(l5q, nh)*2-1)))
		end
		l5qsignal.signal[i] = sig
	end
	# Quantize signal
	if l5qsignal.include_adc
		l5qsignal.signal = round.(l5qsignal.signal.*(2^(nADC-1)-1) ./
                                  maximum(abs.(l5qsignal.signal)))
	end
	return l5qsignal
end
