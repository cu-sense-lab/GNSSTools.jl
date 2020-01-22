include("gnss_tools.jl")
using PyPlot
pygui(true)


"""
    testcourseacquisition(;prn=1, f_d=0., n0=1, t_length=1e-3,
                           fd_range=5000., threads=nthreads(),
                           showplot=false)

Simulates a noisy signal with parameters specified above and
performs course acquisition on it. Prints the course Doppler
and code phase (in samples) estimates. Set `showplot` to `true`
to plot along the time index for estimated Doppler bin.
"""
function testcourseacquisition(;prn=1, f_d=0., n0=1, t_length=1e-3,
                                fd_range=5000., threads=nthreads(),
                                showplot=false)
	# Simulate signal with noise
	type = Val(:l5q)
	f_s = 25e6  # Hz
	f_if = 0.  # Hz
	fd_rate = 0.  # Hz/s
	Tsys = 535.  # K
	CN0 = 45.  # dB*Hz
	ϕ = 0.  # rad
	nADC = 4  # bits
	B = 2.046e7  # Hz
	include_carrier = true
	include_adc = true
	include_noise = true
	include_carrier_amplitude = true
	include_neuman_code = true
	data = definesignal(type::Val{:l5q}, prn, f_s, t_length;
	                    f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
	                    CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
	                    include_carrier=include_carrier,
	                    include_adc=include_adc,
	                    include_noise=include_noise,
	                    include_neuman_code=include_neuman_code,
	                    code_start_idx=n0,
	                    include_carrier_amplitude=include_carrier_amplitude)
	generatesignal!(data)
	# Generate replica signal for cross correlation
	replica = definesignal(type::Val{:l5q}, prn, f_s, t_length;
	                       f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
	                       CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
	                       include_carrier=include_carrier,
	                       include_adc=false,
	                       include_noise=false,
	                       include_neuman_code=include_neuman_code,
	                       code_start_idx=1,
	                       include_carrier_amplitude=false)

	# Perform cross correlation using function
	fd_center = 0.  # Hz
	Δfd = 1/data.t_length  # Hz
	sample_num = data.sample_num
	corr_result = gencorrresult(fd_range, Δfd, sample_num)
	courseacquisition!(corr_result, data, replica, prn;
	                   fd_center=fd_center, fd_range=fd_range,
	                   fd_rate=fd_rate, Δfd=Δfd, threads=threads)
	max_idx = argmax(corr_result)
	fd_est = (fd_center-fd_range) + (max_idx[1]-1)*Δfd
	n0_est = max_idx[2]%Int(f_s*nh_code_length/nh_chipping_rate)
	println("\nPRN $(prn):\nfd = $(fd_est)Hz\nn₀ = $(n0_est) samples")
	if showplot
		figure()
		plot(corr_result[max_idx[1],:], "k-")
		xlabel("n (Samples)")
		ylabel("|replica⋆data|² at Peak Doppler Bin")
		title("PRN $(prn)")
	end
end
