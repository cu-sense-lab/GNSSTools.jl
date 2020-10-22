"""
    Satellite
"""
struct Satellite
    id
    init_orbit
    t
    orbit
    r_ecef
    v_ecef
    a_ecef
end


"""
    Constellation
"""
struct Constellation
    epoch
    plane_num
    satellite_per_plane
    Ω₀
    f₀
    ω
    e
    i
    t_range
    ΔΩ
    Δf
    satellites
end


"""
    define_constellation(plane_num, sat_per_plane, inclination)

Set `ΔΩ` to 30ᵒ if simulating sun-sync constellation, such as Iridium, otherwise,
`ΔΩ` will default to `360/plane_num`.
"""
function define_constellation(a, plane_num, satellite_per_plane, i, t_range;
                              eop=get_iers_eop(:IAU1980), show_plot=false,
                              Ω₀=0., f₀=0., ω=0., e=0., t_start=0., obs_lla=missing,
                              ΔΩ=360/plane_num)
    a = float(a)
    i = i*π/180
    t_range = float.(t_range) ./ (60*60*24) .+ t_start
    ΔΩ = ΔΩ*π/180
    Δf = 2π/satellite_per_plane
    if f₀ == 0.
        f₀ += Δf/satellite_per_plane
    end
    k = 1
    p = Progress(Int(plane_num*satellite_per_plane), 1,
                 "Generating constellation...")
    if show_plot
        fig = figure()
        ax = subplot()
    end
    satellites = Array{Satellite}(undef, plane_num*satellite_per_plane)
    for plane in 1:plane_num
        for sat in 1:satellite_per_plane
            id = k
            Ω = Ω₀ + (plane-1)*ΔΩ
            f = f₀ + (sat-1)*Δf
            init_orbit = init_orbit_propagator(Val(:twobody), 0., a, e, i,
                                               Ω, ω, f)
            orbit, r, v = propagate_to_epoch!(init_orbit, t_range)
            r_ecef = Array{Float64,2}(undef, length(orbit), 3)
            v_ecef = Array{Float64,2}(undef, length(orbit), 3)
            a_ecef = Array{Float64,2}(undef, length(orbit), 3)
            for i in 1:length(orbit)
                orbit_teme = kepler_to_sv(orbit[i])
                orbit_ecef = svECItoECEF(orbit_teme, TEME(), ITRF(), t_start, eop)
                r_ecef[i,:] = [orbit_ecef.r[1], orbit_ecef.r[2], orbit_ecef.r[3]]
                v_ecef[i,:] = [orbit_ecef.v[1], orbit_ecef.v[2], orbit_ecef.v[3]]
                a_ecef[i,:] = [orbit_ecef.a[1], orbit_ecef.a[2], orbit_ecef.a[3]]
            end
            # orbit_teme = [kepler_to_sv(orbit[i]) for i in 1:length(orbit)]
            # orbit_ecef = [svECItoECEF(orbit_teme[i], TEME(), ITRF(), 0., eop)
            #               for i in 1:length(orbit_teme)]
            # r_ecef = [orbit_ecef[1].r[1] orbit_ecef[1].r[2] orbit_ecef[1].r[3]]
            # v_ecef = [orbit_ecef[1].v[1] orbit_ecef[1].v[2] orbit_ecef[1].v[3]]
            # a_ecef = [orbit_ecef[1].a[1] orbit_ecef[1].a[2] orbit_ecef[1].a[3]]
            # for i in 2:length(orbit_ecef)
            #     r_ecef = [r_ecef; [orbit_ecef[i].r[1] orbit_ecef[i].r[2] orbit_ecef[i].r[3]]]
            #     v_ecef = [v_ecef; [orbit_ecef[i].v[1] orbit_ecef[i].v[2] orbit_ecef[i].v[3]]]
            #     a_ecef = [a_ecef; [orbit_ecef[i].a[1] orbit_ecef[i].a[2] orbit_ecef[i].a[3]]]
            # end
            satellites[k] = Satellite(id, init_orbit, t_range, orbit, r_ecef,
                                      v_ecef, a_ecef)
            k += 1
            if show_plot
                scatter3D(r_ecef[1,1], r_ecef[1,2], r_ecef[1,3], s=50)
                plot3D(satellites[k-1].r_ecef[:,1], satellites[k-1].r_ecef[:,2],
                       satellites[k-1].r_ecef[:,3], "k")
            end
            next!(p)
        end
        # if show_plot
        #     plot3D(satellites[k-1].r_ecef[:,1], satellites[k-1].r_ecef[:,2],
        #            satellites[k-1].r_ecef[:,3], "k")
        # end
    end
    if show_plot
        axis("off")
        if ~ismissing(obs_lla)
            obs_ecef = GeodetictoECEF(obs_lla[1], obs_lla[2], obs_lla[3])
            scatter3D(obs_ecef[1], obs_ecef[2], obs_ecef[3], s=300, c="r", marker="*", facecolor="white")
        end
        u = Array(range(0, 2π, length=100))
        v = Array(range(0, 1π, length=100))
        u, v = meshgrid(u, v)
        x = Rₑ.*cos.(u).*sin.(v)
        y = Rₑ.*sin.(u).*sin.(v)
        z = Rₑ.*cos.(v)
        plot_wireframe(x, y, z, color="grey", linestyle=":", rcount=20, ccount=30)
    end
    return Constellation(t_start, plane_num, satellite_per_plane, Ω₀, f₀, ω, e, i,
                         t_range, ΔΩ, Δf, satellites)
end


"""
    doppler_distribution(a, plane_num, satellite_per_plane, i, t_range,
                         obs_lla, sig_freq; eop=get_iers_eop(:IAU1980),
                         Ω₀=0., f₀=0., show_plot=true, ω=0., e=0.,
                         t_start=0., ΔΩ=360/plane_num, min_elevation=5.)
"""
function doppler_distribution(a, plane_num, satellite_per_plane, i, t_range,
                              obs_lla, sig_freq; eop=get_iers_eop(:IAU1980),
                              Ω₀=0., f₀=0., show_plot=true, ω=0., e=0.,
                              t_start=0., ΔΩ=360/plane_num, min_elevation=5.)
    constellation = define_constellation(a, plane_num, satellite_per_plane, i,
                                         t_range; eop=eop, show_plot=show_plot,
                                         Ω₀=Ω₀, f₀=f₀, ω=ω, e=e, t_start=t_start,
                                         obs_lla=missing, ΔΩ=ΔΩ)
    obs_ecef = GeodetictoECEF(obs_lla[1], obs_lla[2], obs_lla[3])
    N = length(t_range)*plane_num*satellite_per_plane
    dopplers = Array{Float64}(undef, N)
    elevations = Array{Float64}(undef, N)
    k = 1
    for satellite in constellation.satellites
        for i in length(satellite.t)
            r = view(satellite.r_ecef, i, :)
            sat_range, azimuth, elevation = calcelevation(r, obs_lla)
            if elevation >= min_elevation
                v = view(satellite.v_ecef, i, :)
                dopplers[k] = calcdoppler(r, v, obs_ecef, sig_freq)
                elevations[k] = elevation
                k += 1
            end
        end
    end
    return (dopplers[1:k-1], elevations[1:k-1])
 end



"""
    plot_satellite_orbit(satellite::Satellite; user_lla=missing)
"""
function plot_satellite_orbit(satellite::Satellite; obs_lla=missing)
    figure()
    ax = subplot()
    # Plot Earth
    u = Array(range(0, 2π, length=100))
    v = Array(range(0, 1π, length=100))
    u, v = meshgrid(u, v)
    x = Rₑ.*cos.(u).*sin.(v)
    y = Rₑ.*sin.(u).*sin.(v)
    z = Rₑ.*cos.(v)
    plot_wireframe(x, y, z, color="grey", linestyle=":", rcount=20, ccount=30)
    # Plot user location
    if ~ismissing(obs_lla)
        obs_ecef = GeodetictoECEF(obs_lla[1], obs_lla[2], obs_lla[3])
        scatter3D(obs_ecef[1], obs_ecef[2], obs_ecef[3], s=300, c="r",
                  marker="*", facecolor="white")
    end
    # Plot satellite
    scatter3D(satellite.r_ecef[1,1], satellite.r_ecef[1,2],
              satellite.r_ecef[1,3], s=50)
    plot3D(satellite.r_ecef[:,1], satellite.r_ecef[:,2],
           satellite.r_ecef[:,3], "k")
    axis("off")
end
