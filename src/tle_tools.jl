"""
    OrbitInfo

Struct containing orbit, start time,
and observer location info.
"""
mutable struct OrbitInfo{T1,T2,T3,T4,T5,T6,T7,T8}
    start_time::T1
    start_time_julian_date::T2
    site_loc_lla::T3
    site_loc_ecef::T4
    tle_file_name::T5
    tle::T6
    orb::T7
    eop::T8
end


"""
    initorbitinfo(source_tle::String, target_tle::String, start_time,
                  site_loc_lla)

Initialize the struct OrbitInfo for multiple satellites. Provide
the file names of the individual TLE files.
"""
function initorbitinfo(source_tle::String, target_tle::String, start_time,
                       site_loc_lla)
    start_time_julian_date = DatetoJD(start_time...)
    site_loc_ecef = GeodetictoECEF(site_loc_lla...)
    tle = [read_tle(source_tle)[1], read_tle(target_tle)[1]]
    orb = [init_orbit_propagator(Val(:sgp4), tle[1], sgp4_gc_wgs84),
           init_orbit_propagator(Val(:sgp4), tle[2], sgp4_gc_wgs84)]
    eop = get_iers_eop(:IAU1980)
    tle_names = [source_tle, target_tle]
    return OrbitInfo(start_time, start_time_julian_date,
                     site_loc_lla, site_loc_ecef,
                     tle_names, tle, orb, eop)
end


"""
    initorbitinfo(source_tle::TLE, target_tle::TLE, start_time,
                  site_loc_lla)

Initialize the struct OrbitInfo for multiple satellites. Provide
the already loaded TLE files as `SatelliteToolbox` `TLE` structs.
"""
function initorbitinfo(source_tle::TLE, target_tle::TLE, start_time,
                       site_loc_lla)
    start_time_julian_date = DatetoJD(start_time...)
    site_loc_ecef = GeodetictoECEF(site_loc_lla...)
    tle = [source_tle, target_tle]
    orb = [init_orbit_propagator(Val(:sgp4), tle[1], sgp4_gc_wgs84),
           init_orbit_propagator(Val(:sgp4), tle[2], sgp4_gc_wgs84)]
    eop = get_iers_eop(:IAU1980)
    tle_names = [source_tle.name, target_tle.name]
    return OrbitInfo(start_time, start_time_julian_date,
                     site_loc_lla, site_loc_ecef,
                     tle_names, tle, orb, eop)
end


"""
    initorbitinfo(source_tle::String, start_time, site_loc_lla)

Initialize the struct OrbitInfo for single satellite. Provide
the file name of the individual TLE file.
"""
function initorbitinfo(source_tle::String, start_time, site_loc_lla)
    start_time_julian_date = DatetoJD(start_time...)
    site_loc_ecef = GeodetictoECEF(site_loc_lla...)
    tle = read_tle(source_tle)[1]
    orb = init_orbit_propagator(Val(:sgp4), tle, sgp4_gc_wgs84)
    eop = get_iers_eop(:IAU1980)
    return OrbitInfo(start_time, start_time_julian_date,
                     site_loc_lla, site_loc_ecef,
                     source_tle, tle, orb, eop)
end


"""
    initorbitinfo(source_tle::TLE, start_time, site_loc_lla)

Initialize the struct OrbitInfo for single satellite. Provide
the already loaded TLE file as a `SatelliteToolbox` `TLE` struct.
"""
function initorbitinfo(source_tle::TLE, start_time, site_loc_lla)
    start_time_julian_date = DatetoJD(start_time...)
    site_loc_ecef = GeodetictoECEF(site_loc_lla...)
    tle = source_tle
    orb = init_orbit_propagator(Val(:sgp4), tle, sgp4_gc_wgs84)
    eop = get_iers_eop(:IAU1980)
    return OrbitInfo(start_time, start_time_julian_date,
                     site_loc_lla, site_loc_ecef,
                     source_tle.name, tle, orb, eop)
end


"""
    calcdoppler(orb::Array, gps_tle, julian_date, eop,
                obs_ecef, sig_freq)

Use `svECItoECEF` to convert satellite state vector
(r and v) from ECI to ECEF. `orb` is a list of
orbit objects initialized using `SatelliteToolbox`.
The first orbit object is the source, `gps_orb`, and
the second is object the signal is reflected off of,
`sat_orb`. Only two objects are used. Anything after
the second index of `orb` is ignored.
"""
function calcdoppler(orb::Array, gps_tle, julian_date, eop,
                     obs_ecef, sig_freq)
    gps_orb = orb[1]
    sat_orb = orb[2]
    # gps_orb = init_orbit_propagator(Val{:sgp4}, gps_tle)
    # sat_orb = init_orbit_propagator(Val{:sgp4}, sat_tle)
    # Propagate orbits to "julian_date"
    gps_orb, rg, vg = propagate_to_epoch!(gps_orb, julian_date)
    sat_orb, rs, vs = propagate_to_epoch!(sat_orb, julian_date)
    # Convert orbits to state vectors
    gps_teme = kepler_to_sv(gps_orb)
    sat_teme = kepler_to_sv(sat_orb)
    # Transform TEME to ECEF frame
    gps_ecef = svECItoECEF(gps_teme, TEME(), ITRF(), julian_date, eop)
    sat_ecef = svECItoECEF(sat_teme, TEME(), ITRF(), julian_date, eop)
    # Get velocities and positions
    rg = gps_ecef.r
    vg = gps_ecef.v
    rs = sat_ecef.r
    vs = sat_ecef.v
    # Calculate observed radial velocities
    # Radial velocity of satellite relative to observer
    v_obs = transpose((rs-obs_ecef)/norm(rs-obs_ecef))*vs
    # Radial velocity of satellite relative to GPS satellite
    v_g2s = transpose((rs-rg)/norm(rs-rg))*vs
    # Radial velocity of GPS satellite relative to satellite
    v_s2g = transpose((rg-rs)/norm(rg-rs))*vg
    # Calculate observed Doppler frequency by observer
    fd_obs = -sig_freq*(v_obs + v_g2s + v_s2g)/c  # Hz
    return fd_obs
end


"""
    calcdoppler(orb::OrbitPropagatorSGP4, gps_tle, julian_date, eop,
                obs_ecef, sig_freq)

Calculates the Doppler frequency for the direct signal case.
"""
function calcdoppler(orb::OrbitPropagator, gps_tle, julian_date, eop,
                     obs_ecef, sig_freq)
    # Propagate orbits to "julian_date"
    orb, rg, vg = propagate_to_epoch!(orb, julian_date)
    # Convert orbits to state vectors
    orb_teme = kepler_to_sv(orb)
    # Transform TEME to ECEF frame
    orb_ecef = svECItoECEF(orb_teme, TEME(), ITRF(), julian_date, eop)
    # Get velocities and positions
    r = orb_ecef.r
    v = orb_ecef.v
    # Direct signal case
    v_obs = transpose((r-obs_ecef)/norm(r-obs_ecef))*v
    # Calculate observed Doppler frequency by observer
    fd_obs = -sig_freq*(v_obs + v_g2s + v_s2g)/c  # Hz
    return v_obs
end
