include("constants.jl")
include("l5q_code_generator.jl")


"""
L5QSignal()

Struct for holding L5Q GNSS signal properties for
signal generation.
"""
mutable struct L5QSignal{T}
	type::T
	prn::T
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

Define properties of locally generated signal
based off its type, PRN, etc. 
"""
function definesignal(type::Val{:l5q}, prn, f_s, t_length;
	                  f_if=0., f_d=0., fd_rate=0.)


