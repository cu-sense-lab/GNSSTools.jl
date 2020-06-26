"""
    dataprocess(data_file, tle_source, tle_target=missing,
                T=1e-3, t_length=4.)

Perform course acquisition and/or tracking over entire data file using
coherent integration `T` with length of data to be processed being defined
by `t_length`.
"""
function dataprocess(data_file, tle_source=missing, tle_target=missing,
                     start_time=missing, site_loc_lla=missing; T=1e-3,
                     t_length=4., output_dir=missing, prn=:all,
                     data_start_t=1e-3)
    # Load data & infer f_s, f_if, and data type from file name
    f_s, f_if, f_center_data, sigtype, data_type = data_info_from_name(data_file)
    start_data_idx = Int(f_s * start_t)+1
    data = loaddata(data_type, data_file, f_s, f_if, t_length;
                    start_data_idx=start_data_idx)
    # Initialize replica struct
    replica = definesignal(type, f_s, replica_t_length)
    # Calculate Doppler bin spacing for course acquisition
    Δfd = 1/replica.t_length  # Hz
    # Allocate space for correlation result
    corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
    # Initialize orbit struct
    if ismissing(tle_source) && ismissing(tle_target)
        use_orbitinfo = false
    elseif ismissing(tle_target)
        orbitinfo = initorbitinfo(tle_source, start_time, site_loc_lla)
        use_orbitinfo = true
    else
        orbitinfo = initorbitinfo(tle_source, tle_target, start_time,
                                  site_loc_lla)
        use_orbitinfo = true
    end
    # Allocate space for t, fd_est, n₀_est arrays
    N = Int(floor(t_length/T))
    t = calctvector(N, T)
    fd_est = Array{Complex{Float64}}(undef, N)
    n₀_est = Array{Complex{Float64}}(undef, N)
    #
end
