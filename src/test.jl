prn = 26
n0 = 1000.
f_d = 800.
fd_rate = 0.
t_length = 1e-3
fd_range = 5000.
threads = nthreads()
M = 4000

# type = Val(:l5q)
# type = Val(:l5i)
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
    RLM = 20
end

if typeof(type) == Val{:l5i}
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
    include_databits = true
    RLM = 10
end

if typeof(type) == Val{:l1ca}
    f_s = 5e6  # Hz
    # f_if = 1.25e6  # Hz
    f_if = 0.  # Hz
    Tsys = 535.  # K
    CN0 = 45.  # dB*Hz
    ϕ = 0.  # rad
    nADC = 4  # bits
    B = 2.046e6  # Hz
    include_carrier = true
    include_adc = true
    include_noise = true
    include_databits = true
    RLM = 10
end

# Load data
file_dir = "/media/Srv3Pool2/by-location/hi/"
# file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
file_name = "hi_e06_20190411_092347_004814_1575.42M_5.0M_USRP4_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L1
file_path = string(file_dir, file_name)
data_type = Val(:sc4)
start_t = 1e-3
data = loaddata(data_type, file_path, f_s, f_if, M*t_length;
                    start_data_idx=Int(f_s * start_t)+1)

# data = definesignal(type, f_s, M*t_length; prn=prn,
#                     f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
#                     CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
#                     include_carrier=include_carrier,
#                     include_adc=include_adc,
#                     include_noise=include_noise,
#                     code_start_idx=n0)
# if (typeof(type) == Val{:l1ca}) | (typeof(type) == Val{:l5i})
#     data.include_databits = include_databits
# end
# generatesignal!(data)

replica = definesignal(type, f_s, 1e-3; prn=prn,
                           f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                           CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
                           include_carrier=include_carrier,
                           include_adc=false,
                           include_noise=false,
                           code_start_idx=1)
replicalong = definesignal(type, f_s, RLM*t_length; prn=prn,
                           f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                           CN0=CN0, ϕ=ϕ, nADC=nADC, B=B,
                           include_carrier=include_carrier,
                           include_adc=false,
                           include_noise=false,
                           code_start_idx=1)
Δfd = 1/replica.t_length  # Hz
# fd_center = round(f_d/Δfd)*Δfd  # Hz
fd_center = 0.  # Hz
corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
courseacquisition!(corr_result, data, replica, prn;
                   fd_center=fd_center, fd_range=fd_range,
                   fd_rate=fd_rate, Δfd=Δfd, threads=threads)
max_idx = argmax(corr_result)
fd_est = (fd_center-fd_range) + (max_idx[1]-1)*Δfd
if typeof(type) == Val{:l5q}
    n0_est = max_idx[2]#%Int(f_s*nh20_code_length/nh_chipping_rate)
elseif typeof(type) == Val{:l5i}
    n0_est = max_idx[2]#%Int(f_s*nh10_code_length/nh_chipping_rate)
else
    n0_est = max_idx[2]
end
results = fineacquisition(data, replicalong, prn, fd_est,
                          n0_est, Val(:fft))
trackresults = trackprn(data, replica, prn, results.ϕ_init,
                        results.fd_est, results.n0_idx_course)
plotresults(trackresults)