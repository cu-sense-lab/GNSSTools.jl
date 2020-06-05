function calcdoppler(sat_tle, gps_tle, julian_date, eop,
                     obs_ecef)
    """
        Use `svECItoECEF` to convert satellite state vector
        (r and v) from ECI to ECEF.
    """
    # Initialize orbits for TLEs
    gps_orb = init_orbit_propagator(Val{:sgp4}, gps_tle)
    sat_orb = init_orbit_propagator(Val{:sgp4}, sat_tle)
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
    fd_obs = -L5_freq*(v_obs + v_g2s + v_s2g)/c  # Hz
    return fd_obs
end
