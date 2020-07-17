"""
    dataprocess(data_file; target_satnum=missing, T=10e-3, t_length=4.,
                output_dir=missing, prns=:all, data_start_t=5e-3,
                fd_range=5000., fd_center=0., Δdays=5, showplot=true,
                site_lla=(40.01, -105.2437, 1655), figsize=missing,
                use_gps_orbit_only=false)

Perform course acquisition and/or tracking over entire data file using
coherent integration `T` with length of data to be processed being defined
by `t_length`.
"""
function dataprocess(data_file; target_satnum=missing, T=1e-3, t_length=4.,
                     output_dir=missing, prns=:all, data_start_t=5e-3,
                     fd_range=5000., fd_center=0., Δdays=5, showplot=true,
                     site_lla=(40.01, -105.2437, 1655), figsize=missing,
                     use_gps_orbit_only=false)
    # Load data & infer f_s, f_if, and data type from file name
    file_info = data_info_from_name(data_file)
    println("-----------------------------------------------------")
    println("                      File Info")
    println("-----------------------------------------------------")
    println("Timestamp (UTC):            $(file_info.timestamp)")
    println("Sampling frequency (Hz):    $(file_info.f_s)")
    println("Signal type:                $(gnsstypes[file_info.sigtype])")
    println("Data center frequency (Hz): $(file_info.sig_freq)")
    println("IF frequency (Hz):          $(file_info.f_if)")
    println("Data type:                  $(gnsstypes[file_info.data_type])")
    println("-----------------------------------------------------")
    start_data_idx = Int(file_info.f_s * data_start_t)+1
    data = loaddata(file_info.data_type, data_file, file_info.f_s,
                    file_info.f_if, T;
                    start_data_idx=start_data_idx)
    # Initialize replica struct
    replica = definesignal(file_info.sigtype, file_info.f_s, T)
    # Calculate Doppler bin spacing for course acquisition
    Δfd = 1/replica.t_length  # Hz
    # Allocate space for correlation result
    corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
    # Get list of PRNs and corresponding NORAD IDs at given time
    if prns == :all
        if (file_info.sigtype == Val(:l5q)) || (file_info.sigtype == Val(:l5i))
            prns = [1, 24, 27, 30, 6, 9, 3, 26, 8, 10, 32]  # all Block IIR-F satellites
        else
            prns = Array(1:32)
        end
    end
    gps_satnums = getGPSSatnums(file_info.timestamp_JD, prns)
    prns = collect(keys(gps_satnums))
    # Initialize orbit struct if target object TLE given
    if ismissing(target_satnum) && ~use_gps_orbit_only
        use_orbitinfo = false
    else
        use_orbitinfo = true
        if use_gps_orbit_only
            satnum_list = collect(values(gps_satnums))
        else
            satnum_list = [collect(values(gps_satnums)); target_satnum]
        end
        # Get TLEs for each PRN and target satnum
        tles = getTLEs(file_info.timestamp_JD,
                       satnum_list,
                       Δdays=Δdays)
        orbitinfo = Dict{Int,OrbitInfo}()
        for prn in prns
            if use_gps_orbit_only
                orbitinfo[prn] = initorbitinfo(tles[gps_satnums[prn]],
                                               file_info.timestamp,
                                               site_lla)
            else
                orbitinfo[prn] = initorbitinfo(tles[gps_satnums[prn]],
                                               tles[target_satnum],
                                               file_info.timestamp,
                                               site_lla)
            end
        end
    end
    # Allocate space for t, fd_est, n₀_est, and SNR_est arrays
    N = Int(floor(t_length/T))
    t = calctvector(N, 1/T)
    results = Dict{Int,Dict{String,Array{Float64}}}()
    for prn in prns
        results[prn] = Dict("n0_est"=>Array{Float64}(undef, N),
                            "fd_est"=>Array{Float64}(undef, N),
                            "SNR_est"=>Array{Float64}(undef, N),
                            "fd_exp"=>Array{Float64}(undef, N))
    end
    # Perform course acquisition on data at every Nᵗʰ T segment for each
    # `prn` in `prns`
    p = Progress(N, 1, "Processing...")
    for n in 1:N
        Δt_JD_low = t[n]/60/60/24  # Days
        Δt_JD_cen = (t[n] + T/2)/60/60/24  # Days
        Δt_JD_high = (t[n] + T)/60/60/24  # Days
        for prn in prns
            # Calculate Doppler at (t + T/2), the middle of the replica signal
            # Assume that the Doppler is constant
            if use_orbitinfo
                if use_gps_orbit_only
                    fd_exp_low = calcdoppler(orbitinfo[prn].orb,
                                             orbitinfo[prn].start_time_julian_date+Δt_JD_low,
                                             orbitinfo[prn].eop,
                                             orbitinfo[prn].site_loc_ecef,
                                             file_info.sig_freq)
                    fd_exp_high = calcdoppler(orbitinfo[prn].orb,
                                              orbitinfo[prn].start_time_julian_date+Δt_JD_high,
                                              orbitinfo[prn].eop,
                                              orbitinfo[prn].site_loc_ecef,
                                              file_info.sig_freq)
                else
                    fd_exp_low = calcdoppler(orbitinfo[prn].orb[1],
                                             orbitinfo[prn].orb[2],
                                             orbitinfo[prn].start_time_julian_date+Δt_JD_low,
                                             orbitinfo[prn].eop,
                                             orbitinfo[prn].site_loc_ecef,
                                             file_info.sig_freq)
                    fd_exp_high = calcdoppler(orbitinfo[prn].orb[1],
                                              orbitinfo[prn].orb[2],
                                              orbitinfo[prn].start_time_julian_date+Δt_JD_high,
                                              orbitinfo[prn].eop,
                                              orbitinfo[prn].site_loc_ecef,
                                              file_info.sig_freq)
                end
                fd_rate = (fd_exp_high - fd_exp_low)/T
            else
                fd_exp_low = fd_center
                fd_rate = 0.
            end
            # Perform course acquisition
            courseacquisition!(corr_result, data, replica, prn;
                               fd_center=fd_exp_low, fd_range=fd_range,
                               fd_rate=fd_rate, Δfd=Δfd, threads=nthreads(),
                               showprogressbar=false)
            # Estimate n₀_est, fd_est, and SNR_est from course acquisition
            # result
            n₀, f_d, snr = course_acq_est(corr_result, fd_exp_low, fd_range, Δfd)
            results[prn]["n0_est"][n] = n₀
            results[prn]["fd_est"][n] = f_d
            results[prn]["fd_exp"][n] = fd_exp_low
            results[prn]["SNR_est"][n] = snr
        end
        # reloaddata!(data, start_data_idx+n*N+1)
        data = loaddata(file_info.data_type, data_file, file_info.f_s,
                        file_info.f_if, T;
                        start_data_idx=start_data_idx+n*N+1)
        next!(p)
    end
    # Save results to HDF5 file
    if ~ismissing(output_dir)
        h5open(output_dir, "w") do file
            for prn in prns
                write(file, "$prn/n0_est", results[prn]["n0_est"])
                write(file, "$prn/fd_est", results[prn]["fd_est"])
                write(file, "$prn/SNR_est", results[prn]["SNR_est"])
                write(file, "$prn/fd_exp", results[prn]["fd_exp"])
            end
        end
    end
    # Plot results for only 1 PRN. More than 1 are not plotted.
    if showplot && (length(prns) == 1)
        if t_length > 60
            t .*= 1/60
            t_label = "Minutes"
        else
            t_label = "Seconds"
        end
        if ~ismissing(figsize)
            fig = figure(figsize=figsize)
        else
            fig = figure()
        end
        ax1 = subplot(3, 1, 1)
        plot(t, results[prns[1]]["n0_est"], "k.")
        xlabel("Time ($t_label)")
        ylabel("n₀ (Samples)")
        ax2 = subplot(3, 1, 2)
        plot(t, results[prns[1]]["fd_est"]./1000, "k.")
        xlabel("Time ($t_label)")
        ylabel("Doppler Frequency (kHz)")
        ax3 = subplot(3, 1, 3)
        plot(t, results[prns[1]]["SNR_est"], "k.")
        xlabel("Time ($t_label)")
        ylabel("Correlation Peak SNR (dB)")
        suptitle("PRN $(prns[1])")
        show()
    end
    return results
end
