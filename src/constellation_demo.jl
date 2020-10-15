"""
    demo(a, plane_num, satellite_per_plane, user_lla; sigtype="l1ca",
         include_carrier=true, include_adc=true, include_noise=true,
         include_databits=true, T="short", showplot=true, prn=26, n0=1000.,
         t_length=1e-3, M=4000, fd_range=5000., dll_b=8., state_num=3,
         dynamickf=true, covMult=1., q_a=100., figsize=missing, CN0=45.,
         plot3d=true, show_acq_plot=true, saveto=missing)
"""
function demo(a, plane_num, satellite_per_plane, user_lla; sigtype="l1ca",
              include_carrier=true, include_adc=true, include_noise=true,
              include_databits=true, T="short", showplot=true, prn=26, n0=1000.,
              T=1e-3, M=4000, fd_range=5000., dll_b=8., state_num=3,
              dynamickf=true, covMult=1., q_a=100., figsize=missing, CN0=45.,
              plot3d=true, show_acq_plot=true, saveto=missing, inclination=90.)
    #
    # Select signal type
    if sigtype == "l5q"
        type = Val(:l5q)
    elseif sigtype == "l5i"
        type = Val(:l5i)
    elseif sigtype == "l1ca"
        type = Val(:l1ca)
    end

    threads = nthreads()

    # L5Q parameters
    if typeof(type) == Val{:l5q}
        f_s = 25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e7  # Hz
        RLM = 20
        file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
    end

    if typeof(type) == Val{:l5i}
        f_s = 25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e7  # Hz
        RLM = 10
        file_name = "hi_e06_20190411_092347_004814_1176.45M_25.0M_USRP8_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L5
    end

    if typeof(type) == Val{:l1ca}
        f_s = 5e6  # Hz
        # f_if = 1.25e6  # Hz
        f_if = 0.  # Hz
        Tsys = 535.  # K
        phi = π/4  # rad
        nADC = 4  # bits
        B = 2.046e6  # Hz
        RLM = 10
        file_name = "hi_e06_20190411_092347_004814_1575.42M_5.0M_USRP4_X300_LB-SJ-10100-SF_Dish-LinZ.sc4"  # L1
    end

    # Simulate data
    print("Generating PRN $(prn) $(sigtype) signal...")
    data = definesignal(type, f_s, M*T; prn=prn,
                        f_if=f_if, f_d=f_d, fd_rate=fd_rate, Tsys=Tsys,
                        CN0=CN0, ϕ=phi, nADC=nADC, B=B,
                        include_carrier=include_carrier,
                        include_adc=include_adc,
                        include_noise=include_noise,
                        code_start_idx=n0)
    if (typeof(type) == Val{:l1ca}) | (typeof(type) == Val{:l5i})
        data.include_databits = include_databits
    end
    println("Done")
    constellation_t = Array(data.t[1]:1:data.t[end])
    doppler_t = Array(data.t[1]-1:t_length:data.t[end]+1)
    constellation = define_constellation(a, plane_num, satellite_per_plane,
                                         inclination*pi/180, constellation_t;
                                         show_plot=true);

    generatesignal!(data; doppler_curve=doppler_curve, doppler_t=doppler_t)
end
