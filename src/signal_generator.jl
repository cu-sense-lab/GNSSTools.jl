"""
    l5q

Array of symbols specifying that a signal
is a an L5 dataless signal.
"""
const l5q = [:l5, :q]


"""
    l5i

Array of symbols specifying that a signal
is a an L5 signal with data bits.
"""
const l5i = [:l5, :i]


"""
    L5QSignal()

Struct for holding L5Q GNSS signal properties for
signal generation.
"""
mutable struct L5QSignal
	type::Val{:l5q}
	prn::Int64
	f_s::Float64
	t_length::Float64
	f_if::Float64
	f_d::Float64
	fd_rate::Float64
	Tsys::Float64
	CN0::Float64
	ϕ::Float64
	nADC::Int64
	B::Float64
	code_start_idx::Float64
	l5q_init_code_phase::Float64
	nh_init_code_phase::Float64
	t::Array{Float64,1}
	data::Array{Complex{Float64},1}
	include_carrier::Bool
	include_adc::Bool
	include_noise::Bool
	f_l5q_d::Float64
	f_l5q_dd::Float64
	f_nh_d::Float64
	f_nh_dd::Float64
	sample_num::Int64
	isreplica::Bool
    noexp::Bool
    chipping_rate::Float64
    sig_freq::Float64
    code_length::Int64
end


"""
    L5ISignal()

Struct for holding L5Q GNSS signal properties for
signal generation.
"""
mutable struct L5ISignal
    type::Val{:l5q}
    prn::Int64
    f_s::Float64
    t_length::Float64
    f_if::Float64
    f_d::Float64
    fd_rate::Float64
    Tsys::Float64
    CN0::Float64
    ϕ::Float64
    nADC::Int64
    B::Float64
    code_start_idx::Float64
    l5q_init_code_phase::Float64
    databit_init_code_phase::Float64
    t::Array{Float64,1}
    data::Array{Complex{Float64},1}
    include_carrier::Bool
    include_adc::Bool
    include_noise::Bool
    include_databits::Bool
    f_l5q_d::Float64
    f_l5q_dd::Float64
    f_db_d::Float64  # databit Doppler
    f_db_dd::Float64
    sample_num::Int64
    isreplica::Bool
    noexp::Bool
    chipping_rate::Float64
    sig_freq::Float64
    code_length::Float64
end


"""
    definesignal(type::Val{:l5q}, prn, f_s, t_length;
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., ϕ=0., nADC=4, B=2.046e7,
                 include_carrier=true, include_adc=true,
                 include_noise=true, code_start_idx=1)

Define properties of locally generated L5Q signal
based off its type, PRN, etc.
"""
function definesignal(type::Val{:l5q}, prn, f_s, t_length;
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=4, B=2.046e7,
                      include_carrier=true, include_adc=true,
                      include_noise=true, code_start_idx=1)
	sample_num = Int(f_s * t_length)
    # Generate time vector
    t = Array{Float64}(undef, sample_num)
    @threads for i in 1:sample_num
        @inbounds t[i] = (i-1)/f_s  # s
    end
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
	data = Array{Complex{Float64}}(undef, sample_num)
	isreplica = false
	return L5QSignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                     Tsys, CN0, ϕ, nADC, B, code_start_idx,
                     l5q_init_code_phase, nh_init_code_phase, t,
                     data, include_carrier, include_adc,
                     include_noise,
                     f_l5q_d, f_l5q_dd,
                     f_nh_d, f_nh_dd, sample_num,
                     isreplica, false, L5_chipping_rate,
                     L5_freq, L5_code_length)
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
                  isreplica=signal.isreplica)

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
                       noexp=signal.noexp)
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
	signal.f_l5q_d = f_l5q_d
	signal.f_l5q_dd = f_l5q_dd
	signal.f_nh_d = f_nh_d
	signal.f_nh_dd = f_nh_dd
	signal.l5q_init_code_phase = l5q_init_code_phase
	signal.nh_init_code_phase = nh_init_code_phase
    signal.code_start_idx = code_start_idx
	signal.isreplica = isreplica
    signal.noexp = noexp
	return signal
end
