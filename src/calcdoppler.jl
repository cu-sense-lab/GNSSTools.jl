"""
    calcdoppler(gps_orb::OrbitPropagator,
                sat_orb::OrbitPropagator,
                julian_date, eop,
                obs_ecef, sig_freq)


Use `svECItoECEF` to convert satellite state vector (r and v) from ECI to ECEF.
`orb` is a list of orbit objects initialized using `SatelliteToolbox`. The first
orbit object is the source, `gps_orb`, and the second is object the signal is
reflected off of, `sat_orb`. Only two objects are used. Anything after the
second index of `orb` is ignored.


Required Arguments:

- `gps_orb::OrbitPropagator`: GPS satellite `OrbitPropagator` struct
- `sat_orb::OrbitPropagator`: second satellite `OrbitPropagator` struct
- `julian_date::Float64`: Julian date
- `eop`: Earth Orientation Parameters data which is required to convert to ECEF
    * refer to `SatelliteToolbox` module for more info
- `obs_ecef::Vector`: receiver ECEF position in meters with format `(x, y, z)`
- `sig_freq::Float64`: signal carrier frequency in Hz


Returns:

- `fd_obs::Float64`: receiver observed Doppler frequency in Hz
"""
function calcdoppler(gps_orb::OrbitPropagator,
                     sat_orb::OrbitPropagator,
                     julian_date, eop,
                     obs_ecef, sig_freq)
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
    fd_obs = -sig_freq*(v_obs + v_g2s + v_s2g)/c  # Hz
    return fd_obs
end


"""
    calcdoppler(r::Vector, v::Vector, obs_ecef, sig_freq)


Calculates the Doppler frequency for the direct signal case.


Required Arguments:

- `rg::Vector`: first object ECEF position vector in meters with format `(x, y, z)`
- `vg::Vector`: first object ECEF velocity vector in meters with format `(v_x, v_y, v_z)`
- `rs::Vector`: second object ECEF position vector in meters with format `(x, y, z)`
- `vs::Vector`: second object ECEF velocity vector in meters with format `(v_x, v_y, v_z)`
- `julian_date::Float64`: Julian date
- `obs_ecef::Vector`: receiver ECEF position in meters with format `(x, y, z)`
- `sig_freq::Float64`: signal carrier frequency in Hz


Returns:

- `fd_obs::Float64`: receiver observed Doppler frequency in Hz
"""
function calcdoppler(rg::Vector, vg::Vector, rs::Vector, vs::Vector,
                     obs_ecef, sig_freq)
    # Calculate observed radial velocities
    # Radial velocity of GPS satellite relative to observer
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
    calcdoppler(orb::OrbitPropagator, julian_date, eop,
                obs_ecef, sig_freq)


Calculates the Doppler frequency for the direct signal case.


Required Arguments:

- `orb::OrbitPropagator`: satellite `OrbitPropagator` struct
- `julian_date::Float64`: Julian date
- `eop`: Earth Orientation Parameters data which is required to convert to ECEF
    * refer to `SatelliteToolbox` module for more info
- `obs_ecef::Vector`: receiver ECEF position in meters with format `(x, y, z)`
- `sig_freq::Float64`: signal carrier frequency in Hz


Returns:

- `fd_obs::Float64`: receiver observed Doppler frequency in Hz
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
    v_obs = transpose(v)*(r-obs_ecef) / norm(r-obs_ecef)
    # Calculate observed Doppler frequency by observer
    fd_obs = -sig_freq*v_obs/c  # Hz
    return fd_obs
end


"""
    calcdoppler(r::Vector, v::Vector, obs_ecef, sig_freq)


Calculates the Doppler frequency for the direct signal case.


Required Arguments:

- `r`: object ECEF position vector in meters with format `(x, y, z)`
- `v`: object ECEF velocity vector in m/s with format `(v_x, v_y, v_z)`
- `obs_ecef`: receiver ECEF position in meters with format `(x, y, z)`
- `sig_freq`: signal carrier frequency in Hz


Returns:

- `fd_obs::Float64`: receiver observed Doppler frequency in Hz
"""
function calcdoppler(r, v, obs_ecef, sig_freq)
    # Direct signal case
    v_obs = transpose(v)*(r-obs_ecef) / norm(r-obs_ecef)
    # Calculate observed Doppler frequency by observer
    fd_obs = -sig_freq*v_obs/c  # Hz
    return fd_obs
end


"""
    calcdopplerrate(r::Vector, v::Vector, a::Vector, obs_ecef, sig_freq)


Calculates the Doppler frequency for the direct signal case.


Required Arguments:

- `r`: object ECEF position vector in meters with format `(x, y, z)`
- `v`: object ECEF velocity vector in m/s with format `(v_x, v_y, v_z)`
- `a`: object ECEF velocity vector in m/s² with format `(a_x, a_y, a_z)`
- `obs_ecef`: receiver ECEF position in meters with format `(x, y, z)`
- `sig_freq`: signal carrier frequency in Hz


Returns:

- `fd_rate_obs::Float64`: receiver observed Doppler frequency in Hz
"""
function calcdopplerrate(r, v, a, obs_ecef, sig_freq)
    # Direct signal case
    a_obs = transpose(a)*(r-obs_ecef) / norm(r-obs_ecef)  + 
            norm(v)^2 / norm(r-obs_ecef) -
            (transpose(v)*(r-obs_ecef) / norm(r-obs_ecef))^2 / norm(r-obs_ecef)
    # Calculate observed Doppler frequency rate by observer
    fd_rate_obs = -sig_freq*a_obs/c  # Hz/s
    return fd_rate_obs
end
