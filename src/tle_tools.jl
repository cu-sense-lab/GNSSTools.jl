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
    lat = deg2rad(site_loc_lla[1])
    lon = deg2rad(site_loc_lla[2])
    height = site_loc_lla[3]
    site_loc_ecef = GeodetictoECEF(lat, lon, height)
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
    lat = deg2rad(site_loc_lla[1])
    lon = deg2rad(site_loc_lla[2])
    height = site_loc_lla[3]
    site_loc_ecef = GeodetictoECEF(lat, lon, height)
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
    lat = deg2rad(site_loc_lla[1])
    lon = deg2rad(site_loc_lla[2])
    height = site_loc_lla[3]
    site_loc_ecef = GeodetictoECEF(lat, lon, height)
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
    lat = deg2rad(site_loc_lla[1])
    lon = deg2rad(site_loc_lla[2])
    height = site_loc_lla[3]
    site_loc_ecef = GeodetictoECEF(lat, lon, height)
    tle = source_tle
    orb = init_orbit_propagator(Val(:sgp4), tle, sgp4_gc_wgs84)
    eop = get_iers_eop(:IAU1980)
    return OrbitInfo(start_time, start_time_julian_date,
                     site_loc_lla, site_loc_ecef,
                     source_tle.name, tle, orb, eop)
end


"""
    SatelliteRAE

Holds information on a given satellite range,
azimuth, and elevation from a observer position.
"""
struct SatelliteRAE{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11}
    name::A1
    sat_tle::A2
    julian_date_range::A3
    obs_lla::A4
    obs_ecef::A5
    Δt::A6
    ts::A7
    sat_range::A8
    sat_azimuth::A9
    sat_elevation::A10
    sat_ecef::A11
end


"""
    calcenumatrix(obs_lla)

Calculate the ECEF to ENU transformation matrix
using the observer's position in LLA.

**NOTE:** Latitudes and longitudes are in radians.
"""
function calcenumatrix(obs_lla)
    lat = obs_lla[1]  # rad
    lon = obs_lla[2]  # rad
    h = obs_lla[3]  # meters
    return [         -sin(lon)           cos(lon)        0;
            -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat);
             cos(lat)*cos(lon)  cos(lat)*sin(lon) sin(lat)]
end


"""
    calcelevation(sat_tle::TLE, julian_date_range, eop, obs_ecef)

Calculates the elevation of a given satellite relative to the
observer for every second between the range specified in
`julian_date_range`.
"""
function calcelevation(sat_tle::TLE, julian_date_range, eop, obs_lla;
                       name="Satellite"; Δt=1/60/60/24)
    obs_ecef = GeodetictoECEF(obs_lla[1], obs_lla[2], obs_lla[3])
    sat_orb = init_orbit_propagator(Val{:sgp4}, sat_tle)
    ts = Array(julian_date_range[1]:Δt:julian_date_range[2])
    # Propagate orbit to ts
    sat_orb, rs, vs = propagate_to_epoch!(sat_orb, ts)
    # Allocate space for storage
    sat_ranges = Array{Float64}(undef, length(ts))
    azs = Array{Float64}(undef, length(ts))
    els = Array{Float64}(undef, length(ts))
    sat_ecefs = Array{Float64}(undef, length(ts), 3)
    # Calculate ENU transformation matrix
    R_ENU = calcenumatrix(obs_lla)
    for i in 1:length(ts)
        # Convert orbits to state vectors
        sat_teme = kepler_to_sv(sat_orb[i])
        # Transform TEME to ECEF frame and extract position
        sat_ecef = svECItoECEF(sat_teme, TEME(), ITRF(), sat_teme.t, eop).r
        # Calculate user-to-sat vector
        user_to_sat = sat_ecef - obs_ecef
        # Caluclate satellite range
        sat_range = norm(user_to_sat)
        # Normalize `user_to_sat`
        user_to_sat_norm = user_to_sat./sat_range
        # Transform normalized user-to-sat vector to ENU
        enu = R_ENU*user_to_sat_norm  # [East, North, South]
        # Calculate satellite azimuth
        az = atan(enu[1], enu[2])
        # Calculate satellite elevation
        el = asin(enu[3]/norm(enu))
        # Save values
        sat_ranges[i] = sat_range
        azs[i] = rad2deg(az)
        els[i] = rad2deg(el)
        sat_ecefs[i,:] = sat_ecef
    end
    return SatelliteRAE(name, sat_tle, julian_date_range,
                        obs_lla, obs_ecef, Δt, ts, sat_ranges,
                        azs, els, sat_ecefs)
end


"""
    calcelevation(sat_tle, julian_date_range, eop, obs_ecef)

Calculates the elevation of a given satellite relative to the
observer for every second between the range specified in
`julian_date_range`.
"""
function calcelevation(satellite::Satellite, obs_lla; name="Satellite")
    obs_ecef = GeodetictoECEF(obs_lla[1], obs_lla[2], obs_lla[3])
    ts = satellite.t
    Δt = ts[2] - ts[1]
    # Allocate space for storage
    sat_ranges = Array{Float64}(undef, length(ts))
    azs = Array{Float64}(undef, length(ts))
    els = Array{Float64}(undef, length(ts))
    # Calculate ENU transformation matrix
    R_ENU = calcenumatrix(obs_lla)
    for i in 1:length(ts)
        sat_ecef = satellite.r_ecef[i]
        # Calculate user-to-sat vector
        user_to_sat = sat_ecef - obs_ecef
        # Caluclate satellite range
        sat_range = norm(user_to_sat)
        # Normalize `user_to_sat`
        user_to_sat_norm = user_to_sat./sat_range
        # Transform normalized user-to-sat vector to ENU
        enu = R_ENU*user_to_sat_norm  # [East, North, South]
        # Calculate satellite azimuth
        az = atan(enu[1], enu[2])
        # Calculate satellite elevation
        el = asin(enu[3]/norm(enu))
        # Save values
        sat_ranges[i] = sat_range
        azs[i] = rad2deg(az)
        els[i] = rad2deg(el)
    end
    return SatelliteRAE(name, sat_tle, julian_date_range,
                        obs_lla, obs_ecef, Δt, ts, sat_ranges,
                        azs, els, satellite.r_ecef)
end


"""
    getGPSSatnums(obs_time_JD, prns)

Get the NORAD ID for each PRN for a given
observationn time in Julian Days. Read
from `GPSData` dictionary. Returns
dictionary with keys being the prns
"""
function getGPSSatnums(obs_time_JD, prns)
    gps_satnums = Dict{Int,Int}()
    for prn in prns
        for sv in collect(keys(GPSData[prn]))
            if GPSData[prn][sv]["active"]
                if obs_time_JD >= GPSData[prn][sv]["start_julian_date"]
                    gps_satnums[prn] = GPSData[prn][sv]["satnum"]
                end
            else
                if ((obs_time_JD > GPSData[prn][sv]["start_julian_date"]) &&
                   (obs_time_JD < GPSData[prn][sv]["end_julian_date"]))
                    gps_satnums[prn] = GPSData[prn][sv]["satnum"]
                end
            end
        end
    end
    return gps_satnums
end


"""
    getTLEs(obs_time_JD, satnums)

Query Space-Track.org for TLEs matching time and satnum
criteria. Parse and determine TLEs for eac satnum that
is closest but before the observation time.
"""
function getTLEs(obs_time_JD, satnums; Δdays=5)
    obs_time_JD_begin = obs_time_JD - Δdays
    obs_time_begin = JDtoDate(obs_time_JD_begin)
    obs_time = JDtoDate(obs_time_JD+1)
    satnum_list = string(satnums[1])
    for i in 2:length(satnums)
        satnum_list = string(satnum_list, ",", string(satnums[i]))
    end
    tle_file = string(homedir(), "/.GNSSTools/tles.tle")
    date_range = string(string(obs_time_begin[1], pad=4), "-",
                        string(obs_time_begin[2], pad=2), "-",
                        string(obs_time_begin[3], pad=2), "--",
                        string(obs_time[1], pad=4), "-",
                        string(obs_time[2], pad=2), "-",
                        string(obs_time[3], pad=2))
    login_url = "https://www.space-track.org/ajaxauth/login"
    base_query = "query=https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/"
    query_tail = "/orderby/NORAD_CAT_ID/format/3le"
    norad_cat_id = "/NORAD_CAT_ID/"
    # Get username from user
    print("Provide Space-Track.org username: ")
    username = readline(stdin)
    secret_pass = Base.getpass("Provide Space-Track.org user password")
    password = read(secret_pass, String)
    run(`curl -o $tle_file $login_url -d identity=$username"""&"""password=$password"""&"""$base_query$date_range$norad_cat_id$satnum_list$query_tail`;
        wait=false);
    Base.shred!(secret_pass)
    tles = read_tle(tle_file)
    filtered_tles = Dict{Int,TLE}()
    for satnum in satnums
        Δts = []
        idxs = []
        for i in 1:length(tles)
            tle = tles[i]
            if tle.sat_num == satnum
                append!(Δts, abs(tle.epoch-obs_time_JD))
                append!(idxs, i)
            end
        end
        if ~isempty(Δts)
            min_Δt_idx = argmin(Δts)
            filtered_tles[satnum] = tles[idxs[min_Δt_idx]]
        else
            @warn "$satnum had no TLE results. Increase Δdays (currently set to $Δdays)."
        end
    end
    return filtered_tles
end
