include("constants.jl")
include("l5q_code_generator.jl")


"""
L5QSignal()

Struct for holding L5Q GNSS signal properties for
signal generation.
"""
mutable struct L5QSignal{T}
	type::T
	prns::T
	f_s::T
	t_length::T
	f_if::T
	f_d::T
	fd_rate::T
	fds::T
	Tsys::T
	CN0::T
	ϕ::T
	nADC::T
	B::T
	code_start_idx::T
	code_phase_offset::T
	t::T
	signal::T
	include_carrier::T
	include_adc::T
	include_noise::T
	l5q_codes::T
	nh_code::T
end


"""
definesignal()

Define properties of locally generated L5Q signal
based off its type, PRN, etc. 
"""
function definesignal(type::Val{:l5q}, prn, f_s, t_length;
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=0, B=2.046e7,
                      include_carrier=true, include_adc=true,
                      include_noise=true, include_neuman_code=true,
                      code_start_idx=1, code_phase_offset=1)
	# Generate time vector
	t = Array(0:1/f_s:t_length-1/f_s)  # s
	# Get codes for signal
	l5q_codes = gen
end


"""
samplecode()


"""
	
