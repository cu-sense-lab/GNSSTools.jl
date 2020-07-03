"""
    dataprocess(data_file, tle_source, tle_target=missing,
                T=1e-3, t_length=4.)

Perform course acquisition and/or tracking over entire data file using
coherent integration `T` with length of data to be processed being defined
by `t_length`.
"""
function dataprocess(data_file; target_satnum=missing, T=10e-3, t_length=4.,
                     output_dir=missing, prns=:all, data_start_t=5e-3,
                     fd_range=5000., fd_center=0., Δdays=3, showplot=true,
                     site_lla=(40.01, -105.2437, 1655))
    # Load data & infer f_s, f_if, and data type from file name
    file_info = data_info_from_name(data_file)
    start_data_idx = Int(file_info.f_s * start_t)+1
    data = loaddata(file_info.data_type, data_file, file_info.f_s,
                    file_info.f_if, t_length;
                    start_data_idx=start_data_idx)
    # Initialize replica struct
    replica = definesignal(file_info.sigtype, file_info.f_s, replica_t_length)
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
    if ismissing(tle_target)
        use_orbitinfo = false
    else
        use_orbitinfo = true
        # Get TLEs for each PRN and target satnum
        tles = getTLEs(file_info.timestamp_JD,
                       [collect(values(gps_satnums)); target_satnum],
                       Δdays=Δdays)
        orbitinfo = Dict{Int,OrbitInfo}()
        for prn in prns
            orbitinfo[prn] = initorbitinfo(tles[gps_satnums[prn]],
                                           tles[target_satnum],
                                           file_info.timestamp_JD,
                                           site_lla)
        end
    end
    # Allocate space for t, fd_est, n₀_est, and SNR_est arrays
    N = Int(floor(t_length/T))
    t = calctvector(N, T)
    results = Dict{Int,Dict{String,Array{Float64}}}()
    for prn in prns
        results[prn] = Dict("n0_est"=>Array{Float64}(undef, N),
                            "fd_est"=>Array{Float64}(undef, N)
                            "SNR_est"=>Array{Float64}(undef, N))
    end
    # Perform course acquisition on data at every Nᵗʰ T segment for each
    # `prn` in `prns`
    for n in 1:N
        Δt_JD = (t[n] + T/2)/60/60/24  # Days
        for prn in prns
            # Calculate Doppler at (t + T/2), the middle of the replica signal
            # Assume that the Doppler is constant
            if use_orbitinfo
                fd_exp = calcdoppler(orbitinfo[prn].orb[1],
                                     orbitinfo[prn].orb[2],
                                     orbitinfo[prn].start_time_julian_date+Δt_JD,
                                     orbitinfo[prn].eop,
                                     orbitinfo[prn].obs_ecef,
                                     file_info.sig_freq)
            else
                fd_exp = fd_center
            end
            # Perform course acquisition
            courseacquisition!(corr_result, data, replica, prn;
                               fd_center=fd_center, fd_range=fd_range,
                               fd_rate=0., Δfd=Δfd, threads=nthreads())
            # Estimate n₀_est, fd_est, and SNR_est from course acquisition
            # result
            n₀, f_d, snr = course_acq_est(corr_result)
            results[prn]["n0_est"][n] = n₀
            results[prn]["fd_est"][n] = f_d
            results[prn]["SNR_est"][n] = snr
        end
    end
    # Save results to HDF5 file
    if ~ismissing(output_dir)
        h5open(output_dir, "w") do file
            for prn in prns
                write(file, "$prn/n0_est", results[prn]["n0_est"])
                write(file, "$prn/fd_est", results[prn]["fd_est"])
                write(file, "$prn/SNR_est", results[prn]["SNR_est"])
            end
        end
    end
    # Plot results for only 1 PRN only
    if showplot && (length(prns) == 1)

    end
end
