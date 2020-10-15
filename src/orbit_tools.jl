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
"""
function define_constellation(a, plane_num, satellite_per_plane, i, t_range;
                              eop=get_iers_eop(:IAU1980), show_plot=false,
                              Ω₀=0., f₀=0., ω=0., e=0., epoch=0., obs_lla=missing)
    a = float(a)
    i = float(i)
    t_range = float.(t_range) ./ (60*60*24) .+ epoch
    ΔΩ = 2π/plane_num
    Δf = 2π/satellite_per_plane
    if f₀ == 0.
        f₀ += Δf/satellite_per_plane
    end
    satellites = []
    k = 1
    p = Progress(Int(plane_num*satellite_per_plane), 1, "Working...")
    if show_plot
        fig = figure()
        ax = subplot()
    end
    for plane in 1:plane_num
        for sat in 1:satellite_per_plane
            id = k
            Ω = Ω₀ + (plane-1)*ΔΩ/2
            f = f₀ + (sat-1)*Δf
            init_orbit = init_orbit_propagator(Val(:twobody), epoch, a, e, i,
                                               Ω, ω, f)
            orbit, r, v = propagate_to_epoch!(init_orbit, t_range)
            orbit_teme = [kepler_to_sv(orbit[i]) for i in 1:length(orbit)]
            orbit_ecef = [svECItoECEF(orbit_teme[i], TEME(), ITRF(), 0., eop)
                          for i in 1:length(orbit_teme)]
            r_ecef = [orbit_ecef[1].r[1] orbit_ecef[1].r[2] orbit_ecef[1].r[3]]
            v_ecef = [orbit_ecef[1].v[1] orbit_ecef[1].v[2] orbit_ecef[1].v[3]]
            a_ecef = [orbit_ecef[1].a[1] orbit_ecef[1].a[2] orbit_ecef[1].a[3]]
            for i in 2:length(orbit_ecef)
                r_ecef = [r_ecef; [orbit_ecef[i].r[1] orbit_ecef[i].r[2] orbit_ecef[i].r[3]]]
                v_ecef = [v_ecef; [orbit_ecef[i].v[1] orbit_ecef[i].v[2] orbit_ecef[i].v[3]]]
                a_ecef = [a_ecef; [orbit_ecef[i].a[1] orbit_ecef[i].a[2] orbit_ecef[i].a[3]]]
            end
            satellite = Satellite(id, init_orbit, t_range, orbit, r_ecef, v_ecef,
                                  a_ecef)
            push!(satellites, satellite)
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
            u = Array(range(0, 2π, length=100))
            v = Array(range(0, 1π, length=100))
            u, v = meshgrid(u, v)
            x = Rₑ.*cos.(u).*sin.(v)
            y = Rₑ.*sin.(u).*sin.(v)
            z = Rₑ.*cos.(v)
            plot_wireframe(x, y, z, color="grey", linestyle=":", rcount=20, ccount=30)
            obs_ecef = GeodetictoECEF(obs_lla[1], obs_lla[2], obs_lla[3])
            scatter3D(obs_ecef[1], obs_ecef[2], obs_ecef[3], s=300, c="r", marker="*", facecolor="white")
        end
    end
    return Constellation(epoch, plane_num, satellite_per_plane, Ω₀, f₀, ω, e, i,
                         t_range, ΔΩ, Δf, satellites)
end
