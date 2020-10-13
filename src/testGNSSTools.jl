"""
    demo(;sigtype="l1ca", include_carrier=true, include_adc=true,
          include_noise=true, include_databits=true, simulatedata=true,
          saveto=missing, T="short", showplot=true, f_d=800.,
          fd_rate=0., prn=26, n0=1000., t_length=1e-3, M=4000,
          fd_range=5000., dll_b=5., state_num=2, dynamickf=true,
          covMult=1., q_a=100., figsize=missing, CN0=45., plot3d=true,
          show_acq_plot=true)

Runs a demo of `GNSSTools` showing major capabilities such as course/fine acquisition
and code/carrier phase and Doppler frequency tracking.
"""
function demo(;sigtype="l1ca", include_carrier=true, include_adc=true,
               include_noise=true, include_databits=true, simulatedata=true,
               saveto=missing, T="short", showplot=true, f_d=800.,
               fd_rate=0., prn=26, n0=1000., t_length=1e-3, M=4000,
               fd_range=5000., dll_b=8., state_num=3, dynamickf=true,
               covMult=1., q_a=100., figsize=missing, CN0=45., plot3d=true,
               show_acq_plot=true, doppler_curve=missing, doppler_t=missing,
               fd_center=missing)
    # Select signal type
    if sigtype == "l5q"
        type = Val(:l5q)
    elseif sigtype == "l5i"
        type = Val(:l5i)
    elseif sigtype == "l1ca"
        type = Val(:l1ca)
    end

    threads = nthreads()

    # L5Q parameters
    if typeof(type) == Val{:l5q}
        f_s = 25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e7  # Hz
        RLM = 20
        file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
    end

    if typeof(type) == Val{:l5i}
        f_s = 25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e7  # Hz
        RLM = 10
        file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
    end

    if typeof(type) == Val{:l1ca}
        f_s = 5e6  # Hz
        # f_if = 1.25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e6  # Hz
        RLM = 10
        file_name = "hi_e06_20190411_092347_004814_1575.42M_5.0M_USRP4_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L1
    end

    if simulatedata
        # Simulate data
        print("Generating PRN $(prn) $(sigtype) signal...")
        data = definesignal(type, f_s, M*t_length; prn=prn,
                            f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                            CN0=CN0, ϕ=phi, nADC=nADC, B=B,
                            include_carrier=include_carrier,
                            include_adc=include_adc,
                            include_noise=include_noise,
                            code_start_idx=n0)
        if (typeof(type) == Val{:l1ca}) | (typeof(type) == Val{:l5i})
            data.include_databits = include_databits
        end
        if ismissing(doppler_t)
            doppler_t = data.t
        end
        generatesignal!(data; doppler_curve=doppler_curve, doppler_t=doppler_t)
        println("Done")
    else
        # Load data
        print("Loading $(sigtype) data...")
        file_dir = "/media/Srv3Pool2/by-location/hi/"
        file_path = string(file_dir, file_name)
        data_type = Val(:sc4)
        start_t = 1e-3
        data = loaddata(data_type, file_path, f_s, f_if, M*t_length;
                        start_data_idx=Int(f_s * start_t)+1)
        println("Done")
    end

    print("Performing course acquisition...")
    if T == "long"
        # Use 20ms coherent integration for L5Q signal and 10ms
        # coherent integration for L5I signal
        if sigtype == "l5i"
            replica_t_length = 10e-3
        elseif sigtype == "l5q"
            replica_t_length = 20e-3
        else
            replica_t_length = 1e-3
        end
    else
        replica_t_length = 1e-3
    end
    # Define 1ms and RLM*1ms signals
    replica = definesignal(type, f_s, replica_t_length)
    replicalong = definesignal(type, f_s, RLM*t_length)
    # Calculate Doppler bin spacing for course acquisition
    Δfd = 1/replica.t_length  # Hz
    if ismissing(fd_center)
        fd_center = round(f_d/Δfd)*Δfd  # Hz
    end
    # fd_center = 0.  # Hz
    # Allocate space for correlation result
    corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
    # Perform course acquisition
    courseacquisition!(corr_result, data, replica, prn;
                       fd_center=fd_center, fd_range=fd_range,
                       fd_rate=fd_rate, Δfd=Δfd, threads=threads)
    n0_est, fd_est, SNR_est = course_acq_est(corr_result, fd_center, fd_range, Δfd)
    println("Done ($(SNR_est)dB)")
    # Perform FFT based fine acquisition
    # Returns structure containing the fine, course,
    # and estimated Doppler frequency
    println("Performing FFT based fine acquisition...")
    results = fineacquisition(data, replicalong, prn, fd_est,
                              n0_est, Val(:fft))
    # Perform code/carrier phase, and Doppler frequency tracking on signal
    # using results from fine acquisition as the intial conditions
    trackresults = trackprn(data, replica, prn, results.phi_init,
                            results.fd_est, results.n0_idx_course,
                            results.P, results.R; DLL_B=dll_b,
                            state_num=state_num, dynamickf=dynamickf,
                            covMult=covMult, qₐ=q_a)
    # Plot results and save if `saveto` is a string
    if showplot
        print("Generating figures...")
        # Course Acquisition Results
        if (sigtype == "l1ca") && show_acq_plot
            if ~plot3d
                fig, ax = subplots(1, 1)
                img = ax.imshow(corr_result, interpolation="none", aspect="auto",
                                extent=[0, size(corr_result)[2], fd_range, -fd_range],
                                vmin=minimum(corr_result), vmax=maximum(corr_result)/4)
                scatter(n0_est, fd_est, color="r", s=400, marker="o",
                        facecolor="none")
            else
                fig = figure()
                ax = Axes3D(fig)
                fd_bins = Array(-fd_range+fd_center:1/replica.t_length:fd_range+fd_center)
                x, y = meshgrid(Array(1:replica.sample_num), fd_bins./1000.)
                surf(x, y, corr_result)
                max_val_idx = argmax(collect(Iterators.flatten(corr_result)))
                scatter3D(collect(Iterators.flatten(x))[max_val_idx],
                          collect(Iterators.flatten(y))[max_val_idx],
                          collect(Iterators.flatten(corr_result))[max_val_idx],
                          c="r")
            end
            xlabel("Samples")
            ylabel("Doppler (Hz)")
            suptitle("PRN $(prn)\nfd course: $(fd_est)Hz; n₀ course: $(n0_est)")
        end
        # Tracking results
        plotresults(trackresults; saveto=saveto, figsize=figsize)
        println("Done")
    end
    return trackresults
end


"""
    courseacquisitiontest(prns="all"; sigtype="l1ca", include_carrier=true,
                          include_adc=true, include_noise=true,
                          include_databits=true, simulatedata=true,
                          T="short", f_d=800., fd_rate=0., n0=1000.,
                          t_length=1e-3, M=4000, fd_range=5000.,
                          sim_prn=26)
"""
function courseacquisitiontest(;prns="all", sigtype="l1ca", include_carrier=true,
                               include_adc=true, include_noise=true,
                               include_databits=true, simulatedata=true,
                               T="short", f_d=800., fd_rate=0., n0=1000.,
                               t_length=1e-3, M=4000, fd_range=5000.,
                               sim_prn=26, CN0=45.)
    # Correlation parameters
    if prns == "all"
        prns = Array(1:32)
    elseif typeof(prns) != Array{Int64,1}
        prns = [prns]
    elseif typeof(prns) == Array{Int64,1}
        # pass
    else
        error("Invalid format for prn. Can be array of PRNs, single prn or `all`.")
    end
    # Select signal type
    if sigtype == "l5q"
        type = Val(:l5q)
    elseif sigtype == "l5i"
        type = Val(:l5i)
    elseif sigtype == "l1ca"
        type = Val(:l1ca)
    end

    threads = nthreads()

    # L5Q parameters
    if typeof(type) == Val{:l5q}
        f_s = 25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e7  # Hz
        RLM = 20
        file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
    end

    if typeof(type) == Val{:l5i}
        f_s = 25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e7  # Hz
        RLM = 10
        file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
    end

    if typeof(type) == Val{:l1ca}
        f_s = 5e6  # Hz
        # f_if = 1.25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e6  # Hz
        RLM = 10
        file_name = "hi_e06_20190411_092347_004814_1575.42M_5.0M_USRP4_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L1
    end

    if simulatedata
        # Simulate data
        print("Generating PRN $(sim_prn) $(sigtype) signal...")
        data = definesignal(type, f_s, M*t_length; prn=sim_prn,
                            f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                            CN0=CN0, ϕ=phi, nADC=nADC, B=B,
                            include_carrier=include_carrier,
                            include_adc=include_adc,
                            include_noise=include_noise,
                            code_start_idx=n0)
        if (typeof(type) == Val{:l1ca}) | (typeof(type) == Val{:l5i})
            data.include_databits = include_databits
        end
        generatesignal!(data)
        println("Done")
    else
        # Load data
        print("Loading $(sigtype) data...")
        file_dir = "/media/Srv3Pool2/by-location/hi/"
        file_path = string(file_dir, file_name)
        data_type = Val(:sc4)
        start_t = 1e-3
        data = loaddata(data_type, file_path, f_s, f_if, M*t_length;
                        start_data_idx=Int(f_s * start_t)+1)
        println("Done")
    end

    print("Performing course acquisition on PRNs:\n$(prns)\n")
    CPUtic();
    for prn in prns
        if T == "long"
            # Use 20ms coherent integration for L5Q signal and 10ms
            # coherent integration for L5I signal
            if sigtype == "l5i"
                replica_t_length = 10e-3
            elseif sigtype == "l5q"
                replica_t_length = 20e-3
            else
                replica_t_length = 1e-3
            end
        else
            replica_t_length = 1e-3
        end
        # Define 1ms and RLM*1ms signals
        replica = definesignal(type, f_s, replica_t_length)
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
        # println("Done")
        # println("n₀: $(n0_est); fd: $(fd_est)Hz")
    end
    CPUtoc()
end
