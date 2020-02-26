prn = 1
n0 = 1000.
f_d = 0.
fd_rate = 0.
t_length = 1e-3
replica_tlength = 1e-3
fd_range = 5000.
threads = nthreads()
M = 4000
# Simulate signal with noise
# type = Val(:l5q)
type = Val(:l1ca)
# L5Q parameters
if typeof(type) == Val{:l5q}
    f_s = 25e6  # Hz
    f_if = 0.  # Hz
    Tsys = 535.  # K
    CN0 = 45.  # dB*Hz
    ϕ = π/4  # rad
    nADC = 4  # bits
    B = 2.046e7  # Hz
    include_carrier = true
    include_adc = true
    include_noise = true
end
if typeof(type) == Val{:l1ca}
    f_s = 5e6  # Hz
    f_if = 1.25e6  # Hz
    Tsys = 535.  # K
    CN0 = 45.  # dB*Hz
    ϕ = 0.  # rad
    nADC = 4  # bits
    B = 2.046e6  # Hz
    include_carrier = false
    include_adc = false
    include_noise = false
    include_databits = false
end
# L1 C/A parameters
data = definesignal(type, f_s, M*t_length; prn=prn,
                    f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                    CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
                    include_carrier=include_carrier,
                    include_adc=include_adc,
                    include_noise=include_noise,
                    code_start_idx=n0)
if typeof(type) == Val{:l1ca}
    data.include_databits = include_databits
end
generatesignal!(data)
replica = definesignal(type, f_s, replica_tlength; prn=prn,
                           f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                           CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
                           include_carrier=include_carrier,
                           include_adc=false,
                           include_noise=false,
                           code_start_idx=1)
replicalong = definesignal(type, f_s, 20*t_length; prn=prn,
                           f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                           CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
                           include_carrier=include_carrier,
                           include_adc=false,
                           include_noise=false,
                           code_start_idx=1)
# Perform cross correlation using function
fd_center = 0.  # Hz
Δfd = 1/replica.t_length  # Hz
corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
courseacquisition!(corr_result, data, replica, prn;
                   fd_center=fd_center, fd_range=fd_range,
                   fd_rate=fd_rate, Δfd=Δfd, threads=threads)
max_idx = argmax(corr_result)
fd_est = (fd_center-fd_range) + (max_idx[1]-1)*Δfd
if typeof(type) == Val{:l5q}
    n0_est = max_idx[2]%Int(f_s*nh_code_length/nh_chipping_rate)
else
    n0_est = max_idx[2]
end
results = fineacquisition(data, replicalong, prn, fd_est,
                          n0_est, Val(:fft))
trackresults = trackprn(data, replica, prn, results.ϕ_init,
                        results.fd_est, results.n0_idx_course)
plotresults(trackresults)