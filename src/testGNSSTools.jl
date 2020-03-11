function demo(;prn=26, file_dir=missing, file_name=missing, sigtype="l5q", threads=nthreads()
               n0=1000., f_d=800., t_length=1e-3, fd_range=5000., M=4000,
               f_if=0., phi=π/4, nADC=4, include_carrier=true, include_adc=true,
               include_noise=true, include_databits=true, start_t=1e-3, B=2.046e7,
               CN0=45., Tsys=535., usesimdata=false, saveto=missing, T=1e-3)
  # Select signal type
  if sigtype == "l5q"
      type = Val(:l5q)
  elseif sigtype == "l5i"
      type = Val(:l5i)
  elseif sigtype == "l1ca"
      type = Val(:l1ca)
  end

  # L5Q parameters
  if typeof(type) == Val{:l5q}
      f_s = 25e6  # Hz
      RLM = 20
      if ismissing(file_name)
        file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
      end
  end

  if typeof(type) == Val{:l5i}
      f_s = 25e6  # Hz
      RLM = 10
      if ismissing(file_name)
        file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
      end
  end

  if typeof(type) == Val{:l1ca}
      f_s = 5e6  # Hz
      RLM = 10
      if ismissing(file_name)
        file_name = "hi_e06_20190411_092347_004814_1575.42M_5.0M_USRP4_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L1
      end
  end

  if usesimdata
      # Simulate data
      data = definesignal(type, f_s, M*t_length; prn=prn,
                          f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                          CN0=CN0, ϕ=phi nADC=nADC, B=B,
                          include_carrier=include_carrier,
                          include_adc=include_adc,
                          include_noise=include_noise,
                          code_start_idx=n0)
      if (typeof(type) == Val{:l1ca}) | (typeof(type) == Val{:l5i})
          data.include_databits = include_databits
      end
      generatesignal!(data)
  else
      # Load data
      if ismissing(file_dir)
        file_dir = "/media/Srv3Pool2/by-location/hi/"
      end
      file_path = string(file_dir, file_name)
      data_type = Val(:sc4)
      data = loaddata(data_type, file_path, f_s, f_if, M*t_length;
                          start_data_idx=Int(f_s * start_t)+1)
  end

  # Define 1ms and RLM*1ms signals
  replica = definesignal(type, f_s, t_length)
  replicalong = definesignal(type, f_s, RLM*t_length)
  # Calculate Doppler bin spacing for course acquisition
  Δfd = 1/replica.t_length  # Hz
  # fd_center = round(f_d/Δfd)*Δfd  # Hz
  fd_center = 0.  # Hz
  # Allocate space for correlation result
  corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
  # Perform course acquisition
  courseacquisition!(corr_result, data, replica, prn;
                     fd_center=fd_center, fd_range=fd_range,
                     fd_rate=fd_rate, Δfd=Δfd, threads=threads)
  # Get peak maximum index location
  max_idx = argmax(corr_result)
  # Calculate course Doppler frequency
  fd_est = (fd_center-fd_range) + (max_idx[1]-1)*Δfd
  # Get course peak index location in time
  n0_est = max_idx[2]
  # Perform FFT based fine acquisition
  # Returns structure containing the fine, course,
  # and estimated Doppler frequency
  results = fineacquisition(data, replicalong, prn, fd_est,
                            n0_est, Val(:fft))
  # Perform code/carrier phase, and Doppler frequency tracking on signal
  # using results from fine acquisition as the intial conditions 
  if ((sigtype == "l5q") | (sigtype == "l5i")) && (T == "long")
      # Use 20ms coherent integration for L5Q signal and 10ms
      # coherent integration for L5I signal
      trackresults = trackprn(data, replicalong, prn, results.phi_init,
                              results.fd_est, results.n0_idx_course)
  else
    # Only use 1ms coherent integration for L1 C/A signal
      trackresults = trackprn(data, replica, prn, results.phi_init,
                              results.fd_est, results.n0_idx_course)
  end
  # Plot results and save if `saveto` is a string
  plotresults(trackresults; saveto=saveto)
end
