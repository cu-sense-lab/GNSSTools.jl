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
