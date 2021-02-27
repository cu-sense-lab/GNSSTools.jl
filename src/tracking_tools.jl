"""
    PLLParms


Struct for storing coefficients for the PLL filter.


Fields:

- `T::Float64`: integration time in seconds
- `damping::Float64`: PLL damping coefficient
- `B::Float64`: Pll bandwidth in Hz
- `wn::Float64`:
- `a0::Float64`:
- `a1::Float64`:
- `a2::Float64`:
- `b0::Float64`:
- `b1::Float64`:
- `b2::Float64`:
- `descriminator::String`: string containing descriminator equation expression
"""
struct PLLParms
    T::Float64
    damping::Float64
    B::Float64
    wn::Float64
    a0::Float64
    a1::Float64
    a2::Float64
    b0::Float64
    b1::Float64
    b2::Float64
    descriminator::String
end


"""
    DLLParms


Struct for storing DLL parameters.


Fields:

- `T::Float64`: integration time in seconds
- `B::Float64`: DLL bandwidth in Hz
- `d::Int64`: correlator sample spacing
	* this is the number of samples between descriminators
	* `NOT` the number of chips
- `descriminator::String`: string containing descriminator equation expression
"""
struct DLLParms
    T::Float64
    B::Float64
    d::Int64
    descriminator::String
end


"""
    TrackResults


A struct containing the parameters for tracking and its results.


Fields:

- `prn::Int64`: PRN number to be tracked
- `signaltype::String`: name of signal type
- `dll_parms::DLLParms`: DLL parameters
- `pll_parms::PLLParms`: PLL parameters
- `M::Int64`: the length of the data over the integration time or the number of
              tracking samples
- `T::Float64`: the tracking integration time in seconds
- `N::Int64`: the number of samples corresponding to a single integration time
- `data_file::String`: name of data file processed
- `data_type::String`: the data type (e.g. `sc4` or `sc8`)
- `data_nADC::Int64`: bit depth of data
- `data_start_idx::Int64`: starting index of data
- `data_t_length::Float64`: length of data in seconds loaded into `GNSSData`
                            struct
- `data_total_t_length::Float64`: total length of data file in seconds
- `data_sample_num::Int64`: number of data samples loaded into `GNSSData` struct
- `data_fs::Float64`: data sampling rate in Hz
- `data_fif::Float64`: data IF frequency in Hz
- `data_start_t::T1`: data start UTC time
- `data_site_lla::T2`: latitude, longitude, and height of receiver location,
                       where data was collected
- `data_init_n0::Float64`: index in data where code starts
- `data_init_code_chip::Float64`: fractional chip number at `t=0s`
- `data_init_phi::Float64`: initial phase in rads
- `data_init_fd::Float64`: initial Doppler frequency in Hz
- `t::Array{Float64,1}`: time vector in seconds with length `M`
- `code_err_meas::Array{Float64,1}`: vector containing the meassured code phase
                                     error
- `code_err_filt::Array{Float64,1}`: vector containing filtered code phase error
- `code_phase_meas::Array{Float64,1}`: vector containing the measured code phase
                                       estimate
- `code_phase_filt::Array{Float64,1}`: vector containing the filtered code phase
                                       estimate
- `n0::Array{Float64,1}`: vector containing the estimated chip per time
- `dphi_meas::Array{Float64,1}`: vector containing the measured phase error
                                 in rads
- `phi::Array{Float64,1}`: vector containing the filtered phase estimate in rads
- `delta_fd::Array{Float64,1}`: vector containing the Doppler frequency error in
                                Hz
- `fds::Array{Float64,1}`: vector containing the filtered Doppler frequency in
                           Hz
- `ZP::Array{Complex{Float64},1}`: vector containing the value of the prompt
                                   correlator at each time step
- `SNR::Array{Float64,1}`: vector containing the SNR of the correlation peak at
                           each time step in dB
- `data_bits::Array{Int64,1}`: extracted databits based off value of prompt
                               correlator
	* if `ZP[i] > 0`, `data_bits[i] = 1`
	* if `ZP[i] < 0`, `data_bits[i] = 0`
	* `NOTE:` there may be more samples corresponding to a single data bit
	          depending on the integration time
- `code_length::Int64`: the length of the primary code sequence


Kalman Filter Parameters:

- `Q::Array{Float64,2}`: process noise matrix `Q`
- `A::Array{Float64,2}`: state transition matrix `A`
- `C::Array{Float64,2}`: measurement matrix `C`
- `P::Array{Float64,2}`: state uncertainties per time step
- `x::Array{Float64,2}`: filted state estimates per time step
- `K::Array{Float64,2}`: KF gain per time step
"""
struct TrackResults{T1,T2}
    prn::Int64
    signaltype::String
    dll_parms::DLLParms
    pll_parms::PLLParms
    M::Int64
    T::Float64
    N::Int64
    data_file::String
    data_type::String
    data_nADC::Int64
    data_start_idx::Int64
    data_t_length::Float64
    data_total_t_length::Float64
    data_sample_num::Int64
    data_fs::Float64
    data_fif::Float64
    data_start_t::T1
    data_site_lla::T2
    data_init_n0::Float64
    data_init_code_chip::Float64
    data_init_phi::Float64
    data_init_fd::Float64
    t::Array{Float64,1}
    code_err_meas::Array{Float64,1}
    code_err_filt::Array{Float64,1}
    code_phase_meas::Array{Float64,1}
    code_phase_filt::Array{Float64,1}
    n0::Array{Float64,1}
    dphi_meas::Array{Float64,1}
    phi::Array{Float64,1}
    delta_fd::Array{Float64,1}
    fds::Array{Float64,1}
    ZP::Array{Complex{Float64},1}
    SNR::Array{Float64,1}
    data_bits::Array{Int64,1}
    code_length::Int64
    # Carrier phase Kalman fiter matrices and results
    Q::Array{Float64,2}
    A::Array{Float64,2}
    C::Array{Float64,2}
    P::Array{Float64,2}
    x::Array{Float64,2}
	K::Array{Float64,2}
end


"""
    definepll(T, B, damping)


Define `PLLParms` struct.


Required Arguments:

- `T`: integration time in seconds
- `B`: PLL bandwidth in Hz
- `ξ`: damping coefficient


Returns:

- `PLLParms` struct
"""
function definepll(T, B, ξ)
    ωₙ = 4*B/(2*ξ + 1/(2*ξ))
    a₀ = ωₙ^2*T^2 + 4*ξ*T*ωₙ + 4
    a₁ = 2*ωₙ^2*T^2 - 8
    a₂ = ωₙ^2*T^2 - 4*ξ*T*ωₙ + 4
    b₀ = ωₙ^2*T^2 + 4*ξ*T*ωₙ
    b₁ = 2*ωₙ^2*T^2
    b₂ = ωₙ^2*T^2 - 4*ξ*T*ωₙ
    descriminator = "ϕ = atan(Imag(ZP)/Real(ZP))"
    return PLLParms(T, ξ, B, ωₙ, a₀, a₁, a₂, b₀, b₁, b₂, descriminator)
end


"""
    definedll(T, B, d)


Define `DLLParms` struct.


Required Arguments:

- `T`: integration time in seconds
- `B`: PLL bandwidth in Hz
- `d`: correlator spacing in number of samples


Returns:

- `DLLParms` struct
"""
function definedll(T, B, d)
    descriminator = "Z4 = 1/d * (abs(ZE) - abs(ZL)) / (abs(ZE) + abs(ZL))"
    return DLLParms(T, B, d, descriminator)
end


"""
    measurephase(ZP)


Calculate the raw phase measurement.


Required Arguments:

- `ZP`: prompt correlator value


Returns:

-
"""
function measurephase(ZP)
    return atan((imag(ZP)/real(ZP)))
end


"""
    filtercarrierphase(pll_parms::PLLParms, ϕ_meas, ϕ_meas_1,
                       ϕ_meas_2, ϕ_filt_1, ϕ_filt_2)


Fiters the raw phase measurement using 2ⁿᵈ order PLL filter.


Required Arguments:

- `pll_parms::PLLParms`: struct containing PLL filter parameters
- `ϕ_meas`: phase measured at `t` step in rads
- `ϕ_meas_1`: phase measured at `t-T` step in rads
- `ϕ_meas_2`: phase measured at `t-2T` step in rads
- `ϕ_filt_1`: filtered phase from `t-T` step in rads
- `ϕ_filt_2`: filtered phase from `t-2T` step in rads


Returns:

- filtered phase from `t` step in rads
"""
function filtercarrierphase(pll_parms::PLLParms, ϕ_meas, ϕ_meas_1,
                            ϕ_meas_2, ϕ_filt_1, ϕ_filt_2)
    a₀ = pll_parms.a0
    a₁ = pll_parms.a1
    a₂ = pll_parms.a2
    b₀ = pll_parms.b0
    b₁ = pll_parms.b1
    b₂ = pll_parms.b2
    return (b₀*ϕ_meas + b₁*ϕ_meas_1 + b₂*ϕ_meas_2 - a₁*ϕ_filt_1 - a₂*ϕ_filt_2)/a₀
end


"""
    Z4(dll_parms::DLLParms, ZE, ZP, ZL)


Calculates the code phase error.


Required Arguments:

- `dll_parms::DLLParms`: struct containing DLL filter parameters
- `ZE`: early correlator value
- `ZP`: prompt correlator value
- `ZL`: late correlator value


Returns:

- unfiltered code phase error
"""
function Z4(dll_parms::DLLParms, ZE, ZP, ZL)
    return 0.5 * (abs(ZE) - abs(ZL)) / (abs(ZE) + abs(ZL))
end


"""
    filtercodephase(dll_parms::DLLParms, current_code_err, last_filt_code_err)

Returns the filtered code phase measusurement.


Required Arguments:

- `dll_parms::DLLParms`: struct containing DLL filter parameters
- `current_code_err`: current code phase error estimated from `Z4` descriminator
- `last_filt_code_err`: filtered code phase error from last time step


Returns:

- filtered code phase error
"""
function filtercodephase(dll_parms::DLLParms, current_code_err, last_filt_code_err)
    return last_filt_code_err+4*dll_parms.T*dll_parms.B*(current_code_err-last_filt_code_err)
end


"""
	shiftandcheck(i, offset, N)


Determines the sample location if shifting the current index `i` by some number,
`offset` in a vector with length `N`. Assumes that the new index can wrap.


Required Arguments:

- `i`: current index
- `offset`: sample offset amount
- `N`: size of vector


Returns:

- `j`: new index location
"""
function shiftandcheck(i, offset, N)
	j = i + offset
	if j > N
		j -= N
	end
	if j < 1
		j  += N
	end
	return j
end


"""
	getcorrelatoroutput!(ZP_array, data, replica, i, N, f_if, f_d,
						 fd_rate, ϕ, d, bin_width=1)


Calculate the early, prompt, and late correlator ouputs. Note that
replica already containts the prompt correlator. Be sure to set
the parameters to `replica` and run `generatesignal!(replica)` before
calling this method.


Required Arguments:

- `ZP_array`: allocated memory for ZP values to use to get SNR of peak
- `data`: vector containin I/Q samples to process
- `replica`: replica generated signal
- `i`: current `iᵗʰ` integration time `T`
- `N`: number of samples in `data` and `replica`
- `f_if`: intermediate frequency (IF) of `data` in Hz
- `f_d`: current Doppler frequency estimate in Hz
- `fd_rate`: current Doppler frequency rate estimate in Hz/s
- `ϕ`: current phase estimate in rads
- `d`: current correlator spacing in sample space


Optional Arguments:

- `bin_width`: half width of frequency range to sum when estimating the SNR
               `(default = 1)`


Returns:

- `ze`: early correlator value
- `zp`: prompt correlator value
- `zl`: late correlator value
- `SNR`: signal-to-noise ratio of frequency domain correlation peak in dB
"""
function getcorrelatoroutput!(ZP_array, data, replica, i, N, f_if, f_d,
	                          fd_rate, ϕ, d, bin_width=1)
    # Initialize correlator results
    ze = 0. + 0im
    zp = 0. + 0im
    zl = 0. + 0im
    datasegment = view(data.data, (i-1)*N+1:i*N)
    ts = view(data.t, 1:N)
	Δf = 1/replica.t_length
    # Perform carrier and phase wipeoff and apply early, prompt, and late correlators
    for j in 1:N
		# Perform carrier wipeoff
		# NOTE: # cis(A) = exp(1im*A)
        @inbounds wipeoff = datasegment[j]*cis(-(2π*(f_if+f_d+0.5*fd_rate*ts[j])*ts[j] + ϕ))
		# Calculate prompt correlator output at specific t
		ZP = conj(replica.data[j]) * wipeoff
        @inbounds zp += ZP
		# Get the index of the ZE and ZL correlators at given t
		zeidx = shiftandcheck(j,  d, N)
		zlidx = shiftandcheck(j, -d, N)
		# Calculate early and late correlator outputs
        @inbounds ze += conj(replica.data[zeidx]) * wipeoff
        @inbounds zl += conj(replica.data[zlidx]) * wipeoff
		# Store ZP result at `j` in `jᵗʰ` index in `ZP_array`
		ZP_array[j] = ZP
    end
	ze = ze/N
    zp = zp/N
    zl = zl/N
	# In-place FFT on ZP_array
	fft!(ZP_array)
	# Find peak of FFT result and sum
	zp_abs_sqrd_max = 0.
	zp_abs_sqrd = 0.
	zp_max_idx = 1
	for j in 1:N
		ZP_abs_sqrd = abs2(ZP_array[j])
		if ZP_abs_sqrd > zp_abs_sqrd_max
			zp_abs_sqrd_max = ZP_abs_sqrd
			zp_max_idx = j
		end
		zp_abs_sqrd += ZP_abs_sqrd
	end
	# Calculate SNR
	max_low = shiftandcheck(zp_max_idx, -bin_width, N)
	max_high = shiftandcheck(zp_max_idx, bin_width, N)
	# PS = zp_abs_sqrd_max
	# PS = abs2(ZP_array[1]) + abs2(ZP_array[N])
	# PN = (zp_abs_sqrd - PS)/(N-2)
	PS = sum(abs2.(ZP_array[max_low:N])) + sum(abs2.(ZP_array[1:max_high]))
	PN = (zp_abs_sqrd - PS)/(N-max_high-(N-max_low))
	SNR = 0.
	try
		SNR = 10*log10(sqrt(PS/PN))
	catch
		SNR = NaN
	end
    return (ze, zp, zl, SNR)
end


"""
	calcA(T, state_num=2)


Calculate the state transition matrix, `A`, which is dependent on the
integration time, `T`. Used in the carrier tracking loop. If `state_num` is set
to `2`, then the KF will track the carrier phase and Doppler. If `state_num` is
set to `3`, then the KF will track the carrier phase, Doppler, and Doppler rate.


Required Arguments:

- `T`: integration time in seconds


Optional Arguments:

- `state_num`: number of states to track `(default = 2)`


Returns:

- either `2x2` or `3x3` state transition matrix
"""
function calcA(T, state_num=2)
	if state_num == 3
		return [1 T T^2/2; 0 1 T; 0 0 1]
	elseif state_num == 2
		return [1 T; 0 1]
	else
		error("Number of states specified must be either 2 or 3.")
	end
end


"""
	calcC(T, state_num=2)


Calculate the measurement matrix, `C`, using the integration time, `T`. Used in
the carrier tracking loop. If `state_num` is set to `2`, then the KF will track
the carrier phase and Doppler. If `state_num` is set to `3`, then the KF will
track the carrier phase, Doppler, and Doppler rate.


Required Arguments:

- `T`: integration time in seconds


Optional Arguments:

- `state_num`: number of states to track `(default = 2)`


Returns:

- either `1x2` or `1x3` measurement matrix
"""
function calcC(T, state_num=2)
	if state_num == 3
		return [1 T/2 T^2/6]
	elseif state_num == 2
		return [1 T/2]
	else
		error("Number of states specified must be either 2 or 3.")
	end
end


"""
	calcQ(T, h₀, h₋₂, qₐ, f_L, state_num=2)


Calculate the process noise covariance matrix, `Q`, which is dependent on the
integration time, `T`, and receiver oscillator h-parameters, `h₀` and `h₋₂`. Can
either specify `state_num` as 2 or 3. `f_L` is the carrier frequency. if
`state_num` is set to `3`, higher order effects from sources such as platform
dynamics, decribed with qₐ, are included. If `state_num` is set to `2`, then
they will not be included.


Required Arguments:

- `T`: integration time in seconds
- `h₀`: white frequency noise
- `h₋₂`: random walk frequency noise
- `qₐ`: platform dynamics in m²/s⁶
- `f_L`: carrier frequency in Hz


Optional Arguments:

- `state_num`: number of states to track `(default = 2)`


Returns:

- either `2x2` or `3x3` process noise matric `Q`
"""
function calcQ(T, h₀, h₋₂, qₐ, f_L, state_num=2)
	qϕ = h₀/2  # oscillator phase PSD
	qω = 2*π^2*h₋₂  # oscillator frequency PSD
	if state_num == 3
		return (2π*f_L)^2 .* [T*qϕ+T^3*qω/3+T^5*qₐ/(20*c^2) T^2*qω/2+T^4*qₐ/(8*c^2) T^3*qₐ/(6*c^2);
		                     T^2*qω/2+T^4*qₐ/(8*c^2) T*qω+T^3*qₐ/(3*c^2) T^2*qₐ/(2*c^2);
							 T^3*qₐ/(6*c^2) T^2*qₐ/(2*c^2) T*qₐ/c^2]
	elseif state_num == 2
		return (2π*f_L)^2 .* [T*qϕ+T^3*qω/3 T^2*qω/2; T^2*qω/2 T*qω]
	else
		error("Number of states specified must be either 2 or 3.")
	end
end


"""
	trackprn(data::GNSSSignal, replica::ReplicaSignal, prn, ϕ_init,
			 fd_init, n0_idx_init, P₀, R; DLL_B=5, PLL_B=15, damping=1.4,
			 fd_rate=0., G=0.2, h₀=1e-21, h₋₂=2e-20, σω=10., qₐ=10.,
			 state_num=2, dynamickf=false, cov_mult=1., q_mult=1.,
			 channel="I")


Perform code and phase tracking on data in `data`.

`replica` decides the signal type. Can pass optional arguments that are minumum
amount to track a given PRN.


Required Arguments:

- `data::GNSSSignal`: either `GNSSData` or `ReplicaSignal` struct
- `replica::ReplicaSignal`: struct to use for replica signal generation
	* `replica.t_length` determines the integration time, `T`, used in tracking
- `prn`: PRN number to track
- `ϕ_init`: initial phase estimate in rads
- `fd_init`: initial Doppler frequency in Hz
- `n0_idx_init`: initial code start index location
- `P₀`: 3x3` diagonal matrix containing the initial uncertainties of the
        carrier phase, Doppler, and Doppler rate
- `R`: initial phase measurement uncertainty


Optional Arguments:

- `DLL_B`: bandwidth of the DLL filter in Hz `(default = 5)`
- `PLL_B`: bandwidth of the PLL filter in Hz `(default = 15)`
- `damping`: PLL filter damping coefficient `(default = 1.4)`
- `fd_rate`: Doppler rate in Hz `(default = 0)`
- `G`: `[DEPRICATED]` carrier phase filter gain
- `h₀`: `(default = 1e-21)`
- `h₋₂`: `(default = 2e-20)`
- `qₐ`: line of site platform dynamics in m²/s⁶ `(default = 1)`
- `state_num`: number of states to track `(default = 3)`
	* if set to `3`, carrier phase, Doppler, and Doppler rate are tracked
	* if set to `2`, only carrier phase and Doppler are tracked
- `dynamickf`: flag to specify if steady state KF gain is used or if KF can
               change from time step to time step `(default = true)`
- `cov_mult`: scalar to inflate the initial covariance matrix `P₀`
              `(default = 1)`
- `q_mult`: scalar to inflate the process noise matrix `Q` `(default = true)`
- `channel`: `[NOT USED]` `(default = "I")`


Returns:

- `TrackResults` struct
"""
function trackprn(data::GNSSSignal, replica::ReplicaSignal, prn, ϕ_init,
	              fd_init, n0_idx_init, P₀, R; DLL_B=5, PLL_B=15, damping=1.4,
				  fd_rate=0., G=0.2, h₀=1e-21, h₋₂=2e-20, qₐ=10.,
				  state_num=2, dynamickf=false, cov_mult=1., q_mult=1.,
				  channel="I")
    # Assign signal specific parameters
    sig_freq = replica.signal_type.sig_freq
	if replica.signal_type.include_I
		chipping_rate = replica.signal_type.I_codes.chipping_rates[1]
    	code_length = replica.signal_type.I_codes.code_lengths[1]
	elseif replica.signal_type.include_Q
		chipping_rate = replica.signal_type.Q_codes.chipping_rates[1]
    	code_length = replica.signal_type.Q_codes.code_lengths[1]
	else
		error("Cannot track since there are no codes defined in signal type.")
	end
    # Compute the spacing between the ZE, ZP, and ZL correlators
	# d = floor(Int, (data.f_s/chipping_rate - 1)/2)
	d = floor(Int, data.f_s/chipping_rate/2)
    # Initialize common variables and initial conditions
    T = replica.t_length
    t_length = data.t_length
    f_s = data.f_s
    f_if = data.f_if
    N = replica.sample_num
    f_d = fd_init
    ϕ = ϕ_init
    f_code_d = chipping_rate*(1. + f_d/sig_freq)
    n0_init = calcinitcodephase(code_length,
                                f_code_d, 0.,
                                f_s, n0_idx_init)
    n0 = n0_init
    M = Int64(floor(data.t_length/replica.t_length))
    t = Array(0:T:M*T-T)
    # Define DLL and PLL parameter structs
    dll_parms = definedll(T, DLL_B, d)
    pll_parms = definepll(T, PLL_B, damping)
    # Allocate array space for tracking results
    code_err_meas = Array{Float64}(undef, M)
    code_err_filt = Array{Float64}(undef, M)
    code_phase_meas = Array{Float64}(undef, M)
    code_phase_filt = Array{Float64}(undef, M)
    n0s = Array{Float64}(undef, M)
    dphi_measured = Array{Float64}(undef, M)
    phi = Array{Float64}(undef, M)
    delta_fd = Array{Float64}(undef, M)
    fds = Array{Float64}(undef, M)
    ZP = Array{Complex{Float64}}(undef, M)
	ZP_array = Array{Complex{Float64}}(undef, N)
    SNR = Array{Float64}(undef, M)
    data_bits = Array{Int64}(undef, M)
    A = calcA(T, state_num)
    C = calcC(T, state_num)
    Q = calcQ(T, h₀, h₋₂, qₐ, sig_freq, state_num) .* q_mult
    P = Array{Float64}(undef, size(A)[1], M)
    x = Array{Float64}(undef, size(A)[1], M)
	K = Array{Float64}(undef, size(A)[1], M)
	FFTW.set_num_threads(1)
	if state_num == 3
    	x⁺ᵢ = [ϕ_init; 2π*(f_if+fd_init); 0.]
		P⁺ᵢ = deepcopy(P₀)
	elseif state_num == 2
		x⁺ᵢ = [ϕ_init; 2π*(f_if+fd_init)]
		P⁺ᵢ = deepcopy(P₀)[1:2,1:2]
	else
		error("Number of states specified must be either 2 or 3.")
	end
	x⁻ᵢ = deepcopy(x⁺ᵢ)
	P⁺ᵢ = cov_mult .* P⁺ᵢ
	Kᵢ = zeros(size(x⁻ᵢ))
	Kfixed = dkalman(A, C, Q, Diagonal(R))
    # p = Progress(M, 1, message)
    # Perform code, carrier phase, and Doppler frequency tracking
    for i in 1:M
		if state_num == 3
			ϕ, ω, ωdot = x⁻ᵢ
		else
			ϕ, ω = x⁻ᵢ
			ωdot = 0.
		end
		f_d = ω/2π - f_if
		fd_rate = ωdot/2π
		f_code_d = chipping_rate*(1. + f_d/sig_freq)
		f_code_dd = chipping_rate*fd_rate/sig_freq
        # Calculate the current code start index
        t₀ = ((code_length-n0)%code_length)/f_code_d
        code_start_idx = t₀*f_s + 1
        # Set signal parameters
        definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, phi=0., f_if=0.,
                      # include_carrier=true,
                      # include_thermal_noise=false,
					  # include_phase_noise=false,
                      # include_adc=false,
                      code_start_idx=code_start_idx,
                      isreplica=true, noexp=true)
        # Generate prompt correlator
        generatesignal!(replica, replica.isreplica)
        # Calculate early, prompt, and late correlator outputs
        ze, zp, zl, snr = getcorrelatoroutput!(ZP_array, data, replica, i, N,
		                                       f_if, f_d, fd_rate, ϕ, d)
        # Estimate code phase error
        n0_err = Z4(dll_parms, ze, zp, zl)  # chips
		# Propogate state uncertaninty
		P⁻ᵢ = A*P⁺ᵢ*A' + Q
        # Measure carrier phase and Doopler frequency errors
		# NOTE: `dϕ_meas` is considered the measurement error, `δy`
        dϕ_meas = measurephase(zp)  # rad
        dfd = dϕ_meas/T  # rad/s
		# Estimate Kalman gain
		Kᵢ = (P⁻ᵢ*C')/(C*P⁻ᵢ*C' + R)
		K[:,i] = Kᵢ
		# Correct state uncertaninty
		P⁺ᵢ = (I - Kᵢ*C)*P⁻ᵢ
        # Correct state
		if dynamickf
			KF_gain = Kᵢ
		else
			KF_gain = Kfixed
		end
        correction = KF_gain.*dϕ_meas
		x⁺ᵢ = x⁻ᵢ + correction
        if i > 1
            # Filter raw code phase error measurement
            n0_err_filtered = filtercodephase(dll_parms, n0_err, code_err_filt[i-1])
        else
            n0_err_filtered = n0_err
        end
        # Save to allocated arrays
        code_err_meas[i] = n0_err
        code_err_filt[i] = n0_err_filtered
        code_phase_meas[i] = (n0 + n0_err)%code_length
        code_phase_filt[i] = (n0 + n0_err_filtered)%code_length
        n0s[i] = n0
        dphi_measured[i] = correction[1]
        phi[i] = x⁺ᵢ[1]
        delta_fd[i] = correction[2]/2π
        fds[i] = x⁺ᵢ[2]/2π - f_if
        ZP[i] = zp
		SNR[i] = snr
        if real(zp) > 0
            data_bits[i] = 1
        else
            data_bits[i] = 0
        end
        # Calculate main code chipping rate at next `i`
        f_code_d = chipping_rate*(1. + (x⁺ᵢ[2]/2π - f_if)/sig_freq)
        # Update code phase with filtered code phase error and propagate to next `i`
		# This essetially makes the DLL rate aided.
        if i > 1
            n0 += (n0_err_filtered + f_code_d*T + 0.5*f_code_dd*T^2)%code_length
        end
		# Propagate x⁺ᵢ to next time step
		x⁻ᵢ = A*x⁺ᵢ
        # next!(p)
    end
	FFTW.set_num_threads(nthreads())
    # Return `TrackResults` struct
    # Specify what to store if `data` is a simulated signal struct
    # or `GNSSData` struct
    if typeof(data) == GNSSData
        file_name = data.file_name
        total_data_length = data.total_data_length
        data_type = data.data_type
        start_data_idx = data.start_data_idx
        site_loc_lla = data.site_loc_lla
        data_start_time = data.data_start_time
    else
        file_name = "N/A"
        total_data_length = data.t_length
        data_type = "N/A"
        start_data_idx = 1
        site_loc_lla = "N/A"
        data_start_time = "N/A"
    end
    return TrackResults(prn, replica.name,
                        dll_parms, pll_parms, M, T, N,
                        file_name, data_type, data.nADC,
                        start_data_idx, data.t_length,
                        total_data_length, data.sample_num, f_s,
                        f_if, data_start_time, site_loc_lla,
                        float(n0_idx_init), n0_init, ϕ_init, fd_init, t,
                        code_err_meas, code_err_filt, code_phase_meas,
                        code_phase_filt, n0s, dphi_measured, phi,
                        delta_fd, fds, ZP, SNR, data_bits, code_length,
                        Q, A, C, P, x, K)
end


"""
	plotresults(results::TrackResults; saveto=missing,
				figsize=missing, doppler_curve=missing,
				doppler_t=missing, CN0=missing)


Plots the tracking results from the `trackprn` method.


Required Arguments:

- `results::TrackResults`: struct containing results from signal tracking


Optional Arguments:

- `saveto`: `String` specifying file name to save figure to `(default = missing)`
- `figsize`: `Tuple` of length 2 used to specify figure size in inches
	* format is `(height, width)`
	* `(default = missing)`
- `doppler_curve`: truth Doppler curve that will be plotted against Doppler
                   estimate `(default = missing)`
- `doppler_t`: time vector that must be included with `doppler_curve`
               `(default = missing)`
- `CN0`: expected carrier-to-noise ratio `C/N₀` value `(default = missing)`


Plots tracking results figure.
"""
function plotresults(results::TrackResults; saveto=missing,
	                 figsize=missing, doppler_curve=missing,
					 doppler_t=missing, CN0=missing)
	if ismissing(figsize)
		figure()
	else
		figure(figsize=figsize)
	end
    matplotlib.gridspec.GridSpec(3,2)
    # Plot code phase errors
    subplot2grid((3,2), (0,0), colspan=1, rowspan=1)
    plot(results.t, results.code_phase_meas.%results.code_length, "k.", label="Measured code phase")
    plot(results.t, results.code_phase_filt.%results.code_length, "b-", label="Filtered code phase")
    xlabel("Time (s)")
    ylabel("Code Phase (chips)")
    title("DLL Tracking")
    legend()
    # Plot filtered and measured phase errors
    subplot2grid((3,2), (0,1), colspan=1, rowspan=1)
    plot(results.t, results.dphi_meas.*180 ./π, "k.", label="Measured Δϕ")
    xlabel("Time (s)")
    ylabel("ϕ (degrees)")
    title("PLL Tracking")
    legend()
	# Doppler frequency estimate
    subplot2grid((3,2), (1,0), colspan=1, rowspan=1)
	if ~ismissing(doppler_curve) && ~ismissing(doppler_t)
		plot(results.t, results.fds, "k.", label="Estimate")
		plot(doppler_t, doppler_curve, "b-", label="Truth")
		legend()
	else
    	plot(results.t, results.fds, "k.")
	end
    xlabel("Time (s)")
    ylabel("Doppler (Hz)")
    title("Doppler Frequency Estimate")
	# SNR estimate
	subplot2grid((3,2), (1,1), colspan=2, rowspan=1)
	if ~ismissing(CN0)
		plot(results.t, results.SNR, "k.", label="Estimate")
		axhline(y=CN0+10*log10(results.T), color="b",
		        label="Truth")
		legend()
	else
    	plot(results.t, results.SNR, "k.")
	end
    xlabel("Time (s)")
    ylabel("SNR (dB)")
    title("SNR Estimate")
    # Plot ZP real and imaginary parts
    subplot2grid((3,2), (2,0), colspan=2, rowspan=1)
    plot(results.t, real(results.ZP), label="real(ZP)")
    plot(results.t, imag(results.ZP), label="imag(ZP)")
    xlabel("Time (s)")
    ylabel("ZP")
    title("Prompt Correlator Output")
    legend()
    suptitle("PRN $(results.prn)\nfd = $(Int(round(results.data_init_fd))) Hz\nn₀ = $(results.data_init_n0) samples")
    # subplots_adjust(hspace=0.4, wspace=0.4)
    subplots_adjust(hspace=0.4, wspace=0.2)
    if ~ismissing(saveto)
        savefig(saveto::String)
    end
end
