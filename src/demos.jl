#-------------------------------------------------------------------------------
#             GNSSTools.jl Signal Simulation and Data Processing Demo
#-------------------------------------------------------------------------------
"""
    demo(;sigtype="l1ca", include_carrier=true, include_adc=true,
          include_noise=true, include_databits=true, simulatedata=true,
          saveto=missing, T="short", showplot=true, f_d=800.,
          fd_rate=0., prn=26, n0=1000., t_length=1e-3, M=4000,
          fd_range=5000., dll_b=2, state_num=3, dynamickf=true,
          cov_mult=1., q_a=1, figsize=missing, CN0=45., plot3d=true,
          show_acq_plot=true, doppler_curve=missing, doppler_t=missing,
          fd_center=missing, sig_freq=missing, signal=missing, q_mult=1,
          print_steps=true, σω=1000)

Runs a demo of `GNSSTools` showing major capabilities such as course/fine acquisition
and code/carrier phase and Doppler frequency tracking on simulated and real
data.
"""
function demo(;sigtype="l1ca", include_carrier=true, include_adc=true,
               include_databits=true, simulatedata=true,
               saveto=missing, T="short", showplot=true, f_d=800.,
               fd_rate=0., prn=26, n0=1000., t_length=30, M=1,
               fd_range=5000., dll_b=5, state_num=3, dynamickf=true,
               cov_mult=1., q_a=1, figsize=missing, CN0=45., plot3d=true,
               show_acq_plot=false, doppler_curve=missing, doppler_t=missing,
               fd_center=0., sig_freq=missing, signal=missing, q_mult=1,
               print_steps=true, σω=1000, channel="I", file_name=missing,
               include_thermal_noise=true, include_phase_noise=true,
               fine_acq_method=:fft, skip_to=0.1, tracking_T=missing,
               save_interval=1, file_location="", file_prefix="SimulatedSignal",
               h_parms=h_parms_tcxo[1])
    if print_steps
        println("Running GNSSTools Signal Simulation and Data Processing Demo")
    end
    # L5Q parameters
    if sigtype == "l5"
        f_s = 25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 8  # bits
        B = 2*L5_chipping_rate  # Hz
        if channel == "I"
            fine_acq_T = 10e-3
        elseif channel == "Q"
            fine_acq_T = 20e-3
        else
            error("Invalid channel specified")
        end
        if ismissing(file_name)
            file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
        end
        f_center = 1176.45e6
        if ismissing(sig_freq)
            sig_freq = L5_freq
        end
        if include_databits
            signal_type = define_l5_code_type(t_length; prns=prn,
                                              sig_freq=sig_freq)
        else
            signal_type = define_l5_code_type(;prns=prn, sig_freq=sig_freq)
        end
    elseif sigtype == "l1ca"
        f_s = 5e6  # Hz
        f_if = 1.25e6  # Hz
        # f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 8  # bits
        B = 2*l1ca_chipping_rate  # Hz
        fine_acq_T = 10e-3
        if ismissing(sig_freq)
            sig_freq = L1_freq
        end
        if ismissing(file_name)
            file_name = "hi_e06_20190411_092347_004814_1575.42M_5.0M_USRP4_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L1
        end
        f_center = 1575.42e6
        if include_databits
            signal_type = define_l1ca_code_type(t_length; prns=prn,
                                                sig_freq=sig_freq)
        else
            signal_type = define_l1ca_code_type(;prns=prn, sig_freq=sig_freq)
        end
    else
        error("Invalid signal type specified.")
    end
    if simulatedata
        # Simulate data
        if ~ismissing(doppler_curve) && ~ismissing(doppler_t)
            zero_idx = findfirst(t -> t == 0, doppler_t)
            f_d = doppler_curve[zero_idx]
            fd_rate = (doppler_curve[zero_idx+1]-doppler_curve[zero_idx]) /
                      (doppler_t[zero_idx+1]-doppler_t[zero_idx])
        end
        if ismissing(signal)
            data = definesignal(signal_type, f_s, save_interval; prn=prn,
                                f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                                CN0=CN0, phi=phi, nADC=nADC,
                                include_carrier=include_carrier,
                                include_adc=include_adc,
                                include_thermal_noise=include_thermal_noise,
                                include_phase_noise=include_phase_noise,
                                code_start_idx=n0,
                                receiver_h_parms=h_parms)
            if print_steps
                p = Progress(floor(Int, t_length/save_interval), 1, "Generating PRN $(prn) $(sigtype) signal...")
            end
            file_name = ""
            operation = "replace"
            if include_thermal_noise
                new_thermal_noise = true
            else
                new_thermal_noise = false
            end
            if include_phase_noise
                new_phase_noise = true
            else
                new_phase_noise = false
            end
            for i in 1:floor(Int, t_length/save_interval)
                generatesignal!(data; doppler_curve=doppler_curve, doppler_t=doppler_t)
                file_name = signal_to_sc_data_file(data, sig_freq, file_location, 
                                                file_prefix;
                                                operation=operation)
                operation = "append"
                phase_noise_init_phase = mean(real.(real(data.phase_noise[end-100:end])))
                definesignal!(data; new_thermal_noise=new_thermal_noise, 
                             new_phase_noise=new_phase_noise,
                             phase_noise_init_phase=phase_noise_init_phase,
                             start_t=i*save_interval)
                if print_steps
                    next!(p)
                end
            end
            data = missing
            GC.gc()
            data = loaddata(:sc8, file_name, f_s, sig_freq-f_if, 
                            sig_freq, t_length)
        else
            data = signal
        end
    else
        # Load data
        if print_steps
            print("Loading $(sigtype) data...")
        end
        file_dir = "/media/Srv3Pool2/by-location/hi/"
        file_path = string(file_dir, file_name)
        data_type = :sc4
        start_t = 1e-3
        data = loaddata(data_type, file_path, f_s, f_center, sig_freq, 
                        t_length; skip_to=skip_to)
        if print_steps
            println("Done")
        end
    end
    if T == "long"
        # Use 20ms coherent integration for L5Q signal and 10ms
        # coherent integration for L5I signal
        if sigtype == "l5"
            if channel == "I"
                acquisition_T = 10e-3
                if ismissing(tracking_T)
                    tracking_T = 10e-3
                end
            elseif channel == "Q"
                acquisition_T = 20e-3
                if ismissing(tracking_T)
                    tracking_T = 20e-3
                end
            else
                error("Invalid signal type specified.")
            end
        else
            acquisition_T = 1e-3
            if ismissing(tracking_T)
                tracking_T = 1e-3
            end
        end
    else
        acquisition_T = 1e-3
        if ismissing(tracking_T)
            tracking_T = 1e-3
        end
    end
    # Process signal
    if print_steps
        print("Processing signal...")
    end
    results = process(data, signal_type, prn, channel;
                      σω=σω, fd_center=fd_center, fd_range=fd_range, M=M,
                      acquisition_T=acquisition_T, fine_acq_T=fine_acq_T,
                      tracking_T=tracking_T, cov_mult=cov_mult,
                      q_a=q_a, q_mult=q_mult, dynamickf=dynamickf, dll_b=dll_b,
                      state_num=state_num, fd_rate=fd_rate, show_plot=false,
                      fine_acq_method=fine_acq_method, return_corrresult=true)
    acqresults, trackresults, corr_result, SNR_est = results
    # Plot results and save if `saveto` is a string
    n0_est = acqresults.n0_idx_course
    fd_est = acqresults.fd_course
    println("Done")
    if showplot
        if print_steps
            print("Generating figures...")
        end
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
        if print_steps
            println("Done")
        end
    end
    return (trackresults, data)
end


"""
    demo(a, plane_num, satellite_per_plane, user_lla=(40.01, -105.2437, 1655);
         sigtype="l1ca", include_carrier=true, include_adc=true, include_noise=true,
         include_databits=true, T="short", showplot=true, prn=26, n0=1000.,
         t_length=1e-3, M=4000, fd_range=5000., dll_b=1., state_num=3,
         dynamickf=true, cov_mult=1., q_a=2., figsize=missing, CN0=45.,
         plot3d=true, show_acq_plot=true, saveto=missing, incl=56.,
         sig_freq=missing, t_start=3/60/24, ΔΩ=360/plane_num, q_mult=1,
         print_steps=true, eop=get_eop(), σω=1000.)

Runs a demo of `GNSSTools` showing major capabilities such as course/fine acquisition
and code/carrier phase and Doppler frequency tracking on simulated data based
off a user defined constellation.

    - `a`: semi-major axis in meters
    - `plane_num`: number of orbital planes in constellation
    - `satellite_per_plane`: number of satellites in each orbital plane
    - `user_lla`: user [latitude, longitude, altitude] in [deg, deg, meters]
                  defaults to (40.01, -105.2437, 1655), or Boulder, CO, USA
"""
function demo(a, plane_num, satellite_per_plane, user_lla=(40.01, -105.2437, 1655);
              sigtype="l1ca", include_carrier=true, include_adc=true, 
              include_noise=true,
              include_databits=true, T="short", showplot=true, prn=26, n0=1000.,
              t_length=60, M=1, fd_range=5000., dll_b=1., state_num=3,
              dynamickf=true, cov_mult=1., q_a=1., figsize=missing, CN0=45.,
              plot3d=true, show_acq_plot=false, saveto=missing, incl=56.,
              sig_freq=missing, t_start=3/60/24, ΔΩ=360/plane_num, q_mult=1,
              print_steps=true, eop=get_eop(), σω=1000., channel="I",
              include_thermal_noise=true, include_phase_noise=true,
              fine_acq_method=:fft, tracking_T=missing, 
              Δf_per_plane=360/satellite_per_plane/2, save_interval=1,
              file_location="", file_prefix="SignalFromLEO",
              h_parms=h_parms_tcxo[1])
    if print_steps
        println("Running GNSSTools Constellation Demo")
    end

    # L5Q parameters
    if sigtype == "l5"
        if ismissing(sig_freq)
            sig_freq = L5_freq
        end
        f_s = 25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 8  # bits
        B = 2.046e7  # Hz
        if channel == "I"
            fine_acq_T = 10e-3
        elseif channel == "Q"
            fine_acq_T = 20e-3
        else
            error("Invalid channel specified")
        end
        if include_databits
            signal_type = define_l5_code_type(t_length; prns=prn,
                                              sig_freq=sig_freq)
        else
            signal_type = define_l5_code_type(;prns=prn, sig_freq=sig_freq)
        end
    elseif sigtype == "l1ca"
        if ismissing(sig_freq)
            sig_freq = L1_freq
        end
        # f_s = 25e6  #
        f_s = 5e6  # Hz
        f_if = 1.25e6  # Hz
        # f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 8  # bits
        B = 2.046e6  # Hz
        fine_acq_T = 10e-3
        if include_databits
            signal_type = define_l1ca_code_type(t_length; prns=prn,
                                                sig_freq=sig_freq)
        else
            signal_type = define_l1ca_code_type(;prns=prn, sig_freq=sig_freq)
        end
    else
        error("Invalid signal type specified.")
    end

    # Simulate data
    if print_steps
        print("Setting signal parameters and generating constellation...")
    end
    data = definesignal(signal_type, f_s, save_interval; prn=prn,
                        f_if=f_if, f_d=0., fd_rate=0., Tsys=Tsys,
                        CN0=CN0, phi=phi, nADC=nADC,
                        include_carrier=include_carrier,
                        include_adc=include_adc,
                        include_thermal_noise=include_thermal_noise,
                        include_phase_noise=include_phase_noise,
                        code_start_idx=n0,
                        receiver_h_parms=h_parms)
    constellation_t = Array(data.start_t:1:t_length)
    if t_length >= 1
        t_step = 0.25
    else
        t_step = 0.001
    end
    doppler_t = Array(data.start_t:t_step:t_length+t_step)
    constellation = define_constellation(a, plane_num, satellite_per_plane,
                                         incl, constellation_t;
                                         show_plot=false, obs_lla=user_lla,
                                         eop=eop, t_start=t_start, ΔΩ=ΔΩ,
                                         Δf_per_plane=Δf_per_plane)
    elevations = []
    max_elevations = []
    for satellite in constellation.satellites
        elevation = calcelevation(satellite, user_lla;
                                  name=string(satellite.id))
        push!(max_elevations, maximum(elevation.sat_elevation))
        push!(elevations, elevation)
    end
    max_idx = argmax(max_elevations)
    user_ecef = GeodetictoECEF(deg2rad(user_lla[1]), deg2rad(user_lla[2]),
                               user_lla[3])
    doppler_curve = zeros(length(doppler_t))
    for i in 1:length(doppler_curve)
        julian_date = doppler_t[i]/60/60/24 .+ t_start
        doppler_curve[i] = calcdoppler(constellation.satellites[max_idx].init_orbit,
                                       julian_date, eop, user_ecef, sig_freq)
    end
    sat_rae = calcelevation(constellation.satellites[max_idx], user_lla)
    sat_elevations = sat_rae.sat_elevation
    sat_azimuth = sat_rae.sat_azimuth
    fd_rate = (doppler_curve[2]-doppler_curve[1])/0.001
    definesignal!(data; f_d=doppler_curve[1], fd_rate=fd_rate)
    if print_steps
        println("Done")
    end
    if print_steps
        p = Progress(floor(Int, t_length/save_interval), 1, "Generating PRN $(prn) $(sigtype) signal...")
    end
    file_name = ""
    operation = "replace"
    if include_thermal_noise
        new_thermal_noise = true
    else
        new_thermal_noise = false
    end
    if include_phase_noise
        new_phase_noise = true
    else
        new_phase_noise = false
    end
    for i in 1:floor(Int, t_length/save_interval)
        generatesignal!(data; doppler_curve=doppler_curve, doppler_t=doppler_t)
        file_name = signal_to_sc_data_file(data, sig_freq, file_location, 
                                           file_prefix;
                                           operation=operation)
        operation = "append"
        phase_noise_init_phase = mean(real.(real(data.phase_noise[end-100:end])))
        definesignal!(data; new_thermal_noise=new_thermal_noise, 
                      new_phase_noise=new_phase_noise,
                      phase_noise_init_phase=phase_noise_init_phase,
                      start_t=i*save_interval)
        if print_steps
            next!(p)
        end
    end
    data = missing
    GC.gc()
    data = loaddata(:sc8, file_name, f_s, sig_freq-f_if, sig_freq, t_length)
    if T == "long"
        # Use 20ms coherent integration for L5Q signal and 10ms
        # coherent integration for L5I signal
        if sigtype == "l5"
            if channel == "I"
                acquisition_T = 10e-3
                if ismissing(tracking_T)
                    tracking_T = 10e-3
                end
            elseif channel == "Q"
                acquisition_T = 20e-3
                if ismissing(tracking_T)
                    tracking_T = 20e-3
                end
            else
                error("Invalid signal type specified.")
            end
        else
            acquisition_T = 1e-3
            if ismissing(tracking_T)
                tracking_T = 1e-3
            end
        end
    else
        acquisition_T = 1e-3
        if ismissing(tracking_T)
            tracking_T = 1e-3
        end
    end
    # Process signal
    if print_steps
        print("Processing signal...")
    end
    Δfd = 1/acquisition_T  # Hz
    fd_center = round(doppler_curve[1]/Δfd)*Δfd  # Hz
    results = process(data, signal_type, prn, channel;
                      σω=σω, fd_center=fd_center, 
                      fd_range=fd_range, M=M,
                      acquisition_T=acquisition_T, 
                      fine_acq_T=fine_acq_T,
                      tracking_T=tracking_T,
                      cov_mult=cov_mult,
                      q_a=q_a, q_mult=q_mult, dynamickf=dynamickf, dll_b=dll_b,
                      state_num=state_num, fd_rate=fd_rate, show_plot=false,
                      fine_acq_method=fine_acq_method, return_corrresult=true)
    acqresults, trackresults, corr_result, SNR_est = results
    # Plot results and save if `saveto` is a string
    n0_est = acqresults.n0_idx_course
    fd_est = acqresults.fd_course
    println("Done")
    # Plot results and save if `saveto` is a string
    if showplot
        if print_steps
            print("Generating figures...")
        end
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
            ylabel("Doppler (kHz)")
            suptitle("PRN $(prn)\nfd course: $(fd_est)Hz; n₀ course: $(n0_est)")
        end
        # Tracking results
        A = sqrt(2*GNSSTools.k*Tsys)*10^(CN0/20)
        σ = sqrt(GNSSTools.k*B*Tsys)
        truth_SNR = 10log10(f_s * tracking_T * A^2 / σ^2)
        plotresults(trackresults; saveto=saveto, figsize=figsize,
                    doppler_curve=doppler_curve, doppler_t=doppler_t,
                    truth_SNR=truth_SNR)
        if print_steps
            println("Done")
        end
    end
    return (trackresults, data, doppler_curve, doppler_t, sat_elevations, 
            sat_azimuth)
end
