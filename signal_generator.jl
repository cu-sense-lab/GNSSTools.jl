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
                         A27,A28,A29}
	type::A1
	prn::A2
	f_s::A3
	t_length::A4
	f_if::A5
	f_d::A6
	fd_rate::A7
	fds::A8
	Tsys::A9
	CN0::A10
	ϕ::A11
	nADC::A12
	B::A13
	code_start_idx::A14
	l5q_init_code_phase::A15
	nh_init_code_phase::A16
	t::A17
	signal::A18
	include_carrier::A19
	include_adc::A20
	include_noise::A21
	include_neuman_code::A22
	l5q_codes::A23
	nh_code::A24
	f_l5q_d::A25
	f_l5q_dd::A26
	f_nh_d::A27
	f_nh_dd::A28
	sample_num::A29
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
                      CN0=45., ϕ=0., nADC=8, B=2.046e7,
                      include_carrier=true, include_adc=true,
                      include_noise=true, include_neuman_code=true,
                      code_start_idx=1, fds=missing)
	# Generate time vector
	t = Array(0:1/f_s:t_length-1/f_s)  # s
	sample_num = Int(f_s * t_length)
	# Get codes for signal
	l5q_codes = gen_l5q_codes(XB_intial_conditions)
	# Neuman sequence
	nh_code = deepcopy(nh20)
	# Calculate code chipping rates with Doppler applied
	# L5Q
	f_l5q_d = L5_chipping_rate*(1. + f_d/L5_freq)
	f_l5q_dd = L5_chipping_rate*fd_rate/L5_freq
	# Neuman sequence
	f_nh_d = nh_chipping_rate*(1. + f_d/L5_freq)
	f_nh_dd = nh_chipping_rate*fd_rate/L5_freq
	# Calculate the L5Q and nh code phase offsets
	t₀ = code_start_idx/f_s
	l5q_init_code_phase = L5_code_length - f_l5q_d*t₀ - 0.5*f_l5q_dd*t₀^2
	nh_init_code_phase = nh_code_length - f_nh_d*t₀ - 0.5*f_nh_dd*t₀^2
	# Allocate space for signal
	signal = Array{Complex{Float64}}(undef, sample_num)
	return L5QSignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                     fds, Tsys, CN0, ϕ, nADC, B, code_start_idx,
                     l5q_init_code_phase, nh_init_code_phase, t,
                     signal, include_carrier, include_adc,
                     include_noise, include_neuman_code,
                     l5q_codes, nh_code, f_l5q_d, f_l5q_dd,
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
function definesignal!(l5qsignal::L5QSignal, type::Val{:l5q}=l5qsignal.type;
                       prn=l5qsignal.prn, f_d=l5qsignal.f_d,
                       f_if=l5qsignal.f_if, fd_rate=l5qsignal.fd_rate,
                       Tsys=l5qsignal.Tsys, CN0=l5qsignal.CN0,
                       ϕ=l5qsignal.ϕ, nADC=l5qsignal.nADC,
                       B=l5qsignal.B,
                       include_carrier=l5qsignal.include_carrier,
                       include_adc=l5qsignal.include_adc,
                       include_noise=l5qsignal.include_noise,
                       include_neuman_code=l5qsignal.include_neuman_code,
                       code_start_idx=l5qsignal.code_start_idx,
                       fds=l5qsignal.fds)
	## Calculate code chipping rates with Doppler applied
	# L5Q
	f_l5q_d = L5_chipping_rate*(1. + f_d/L5_freq)
	f_l5q_dd = L5_chipping_rate*fd_rate/L5_freq
	# Neuman sequence
	f_nh_d = nh_chipping_rate*(1. + f_d/L5_freq)
	f_nh_dd = nh_chipping_rate*fd_rate/L5_freq
	# Calculate the L5Q and nh code phase offsets
	t₀ = code_start_idx/l5qsignal.f_s
	l5q_init_code_phase = L5_code_length - f_l5q_d*t₀ - 0.5*f_l5q_dd*t₀^2
	nh_init_code_phase = nh_code_length - f_nh_d*t₀ - 0.5*f_nh_dd*t₀^2
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
	l5qsignal.fds = fds
	l5qsignal.f_l5q_d = f_l5q_d
	l5qsignal.f_l5q_dd = f_l5q_dd
	l5qsignal.f_nh_d = f_nh_d
	l5qsignal.f_nh_dd = f_nh_dd
	l5qsignal.l5q_init_code_phase = l5q_init_code_phase
	l5qsignal.nh_init_code_phase = nh_init_code_phase
	return l5qsignal
end



"""
generatesignal()

Generates local GNSS signal using paramters defined in a
L5QSignal struct.
"""

