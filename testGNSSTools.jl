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
                                showplot=false, fd_rate=0.)
	# Simulate signal with noise
	type = Val(:l5q)
	f_s = 25e6  # Hz
	f_if = 0.  # Hz
	Tsys = 535.  # K
	CN0 = 45.  # dB*Hz
	ϕ = 0.  # rad
	nADC = 4  # bits
	B = 2.046e7  # Hz
	include_carrier = true
	include_adc = true
	include_noise = true
	data = definesignal(type::Val{:l5q}, prn, f_s, t_length;
	                    f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
	                    CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
	                    include_carrier=include_carrier,
	                    include_adc=include_adc,
	                    include_noise=include_noise,
	                    code_start_idx=n0)
	generatesignal!(data)
	# Generate replica signal for cross correlation
	replica = definesignal(type::Val{:l5q}, prn, f_s, t_length;
	                       f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
	                       CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
	                       include_carrier=include_carrier,
	                       include_adc=false,
	                       include_noise=false,
	                       code_start_idx=1)

	# Perform cross correlation using function
	fd_center = f_d  # Hz
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


"""
    testcourseacquisition(;prn=1, f_d=0., n0=1, t_length=1e-3,
                           fd_range=5000., threads=nthreads(),
                           showplot=false)

Simulates a noisy signal with parameters specified above and
performs course acquisition on it. Prints the course Doppler
and code phase (in samples) estimates. Set `showplot` to `true`
to plot along the time index for estimated Doppler bin.
"""
function testcourseacquisitiondata(;prns=26, t_length=1e-3,
                                    threads=nthreads(),
                                    showplot=false,
                                    fd_center=0.,
                                    fd_range=5000.,
                                    file="test",
                                    start_t=0.)
	# Load data
	f_s = 25e6  # Hz
	f_if = 0.  # Hz
	if file == "test"
		data_type = Val(:sc4)
		file_dir = "/media/Srv3pool2/by-location/hi/"
		file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP5_X300_LB-SJ-10100-SF_Dish-LinW.sc4"
	elseif file == "iss"
		data_type = Val(:sc8)
		file_dir = "/media/share/taylors6/data/"
		file_name = "20191023_195000_n200_cu-dish_1176.45_25.sc8"
	end
	file_path = string(file_dir, file_name)
	data = loaddata(data_type, file_path, f_s, f_if, t_length;
                    start_data_idx=Int(f_s * start_t)+1)
	# Generate replica signal for cross correlation
	replica = definesignal(Val(:l5q), 1, f_s, t_length)
	# Perform cross correlation using function
	fd_rate = 0.  # Hz
	Δfd = 1/t_length  # Hz
	sample_num = data.sample_num
	corr_result = gencorrresult(fd_range, Δfd, sample_num)
	# Set PRNs to course acquire
	if prns == "all"
		prns = Array(1:32)
	elseif typeof(prns) != Array{Int64,1}
		prns = [prns]
	elseif typeof(prns) == Array{Int64,1}
		# pass
	else
		error("Invalid format for prn. Can be array of PRNs, single prn or `all`.")
	end
	# Begin course acquisition
	for prn in prns
		courseacquisition!(corr_result, data, replica, prn;
		                   fd_center=fd_center, fd_range=fd_range,
		                   fd_rate=fd_rate, Δfd=Δfd, threads=threads)
		max_idx = argmax(corr_result)
		fd_est = (fd_center-fd_range) + (max_idx[1]-1)*Δfd
		n0_est = max_idx[2]#%Int(f_s*nh_code_length/nh_chipping_rate)
		snr_est = calcsnr(corr_result[max_idx[1],:])
		println("\nPRN $(prn):")
		println("fd = $(fd_est)Hz")
		println("n₀ = $(n0_est) samples")
		println("SNR = $(snr_est)dB")
		if showplot
			figure()
			plot(data.t.*1000, corr_result[max_idx[1],:], "k-")
			xlabel("Time (ms)")
			ylabel("|replica⋆data|² at Peak Doppler Bin")
			title("PRN $(prn)")
		end
	end
end


