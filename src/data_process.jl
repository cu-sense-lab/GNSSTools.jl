"""
    dataprocess(data_file, tle_source, tle_target=missing,
                T=1e-3, t_length=4.)

Perform course acquisition and/or tracking over entire data file using
coherent integration `T` with length of data to be processed being defined
by `t_length`.
"""
function dataprocess(data_file; target_satnum=missing, T=10e-3, t_length=4.,
                     output_dir=missing, prns=:all, data_start_t=1e-3,
                     fd_range=5000., Δdays=3, showplot=true,
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
    fd_est = Array{Complex{Float64}}(undef, N)
    n₀_est = Array{Complex{Float64}}(undef, N)
    SNR_est = Array{Complex{Float64}}(undef, N)
    #
end
