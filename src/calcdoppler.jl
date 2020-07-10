"""
    calcdoppler(gps_orb::OrbitPropagator,
                sat_orb::OrbitPropagator,
                julian_date, eop,
                obs_ecef, sig_freq)

Use `svECItoECEF` to convert satellite state vector
(r and v) from ECI to ECEF. `orb` is a list of
orbit objects initialized using `SatelliteToolbox`.
The first orbit object is the source, `gps_orb`, and
the second is object the signal is reflected off of,
`sat_orb`. Only two objects are used. Anything after
the second index of `orb` is ignored.
"""
function calcdoppler(gps_orb::OrbitPropagator,
                     sat_orb::OrbitPropagator,
                     julian_date, eop,
                     obs_ecef, sig_freq)
    # gps_orb = init_orbit_propagator(Val{:sgp4}, gps_tle)
    # sat_orb = init_orbit_propagator(Val{:sgp4}, sat_tle)
    # Propagate orbits to "julian_date"
    gps_orbit, rg, vg = propagate_to_epoch!(gps_orb, julian_date)
    sat_orbit, rs, vs = propagate_to_epoch!(sat_orb, julian_date)
    # Convert orbits to state vectors
    gps_teme = kepler_to_sv(gps_orbit)
    sat_teme = kepler_to_sv(sat_orbit)
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
    fd_obs = sig_freq*(v_obs + v_g2s + v_s2g)/c  # Hz
    return -fd_obs
end


"""
    calcdoppler(orb::OrbitPropagatorSGP4, julian_date, eop,
                obs_ecef, sig_freq)

Calculates the Doppler frequency for the direct signal case.
"""
function calcdoppler(orb::OrbitPropagator, julian_date, eop,
                     obs_ecef, sig_freq)
    # Propagate orbits to "julian_date"
    orbit, rg, vg = propagate_to_epoch!(orb, julian_date)
    # Convert orbits to state vectors
    orb_teme = kepler_to_sv(orbit)
    # Transform TEME to ECEF frame
    orb_ecef = svECItoECEF(orb_teme, TEME(), ITRF(), julian_date, eop)
    # Get velocities and positions
    r = orb_ecef.r
    v = orb_ecef.v
    # Direct signal case
    v_obs = transpose((r-obs_ecef)/norm(r-obs_ecef))*v
    # Calculate observed Doppler frequency by observer
    fd_obs = sig_freq*v_obs/c  # Hz
    return -fd_obs
end
