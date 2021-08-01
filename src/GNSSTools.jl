__precompile__()

module GNSSTools

    include("greeting.jl")
    export logo, greeting
    # printstyled(logo, color=:default, bold=false)
    # println(greeting)

    # using Reexport

    using Random
    using FFTW
    using Base.Threads
    # Set the number of threads for FFTW to use. It is set to the number of
    # threads that Julia is set to use. User should launch Julia with the number
    # of threads they wish to use.
    FFTW.set_num_threads(nthreads())
    using ProgressMeter
    using Statistics
    using LinearAlgebra
    using CPUTime
    using SatelliteToolbox
    using Interpolations
    using QuadGK
    using HDF5
    using StatsBase
    using SpecialFunctions
    using Optim
    using PyPlot

    pygui(true)

    include("constants.jl")
    include("general_tools.jl")
    include("signal_types.jl")
    include("generate_phase_noise.jl")
    include("calc_code_val.jl")
    include("definesignal.jl")
    include("generatesignal.jl")
    include("codes/l1ca_code_generator.jl")
    include("codes/l5_code_generator.jl")
    include("gnss_data_tools.jl")
    include("course_acquisition_tools.jl")
    include("fine_acquisition_tools.jl")
    include("rinex_tools.jl")
    include("tracking_tools.jl")
    include("dkalman.jl")
    include("calcdoppler.jl")
    include("orbit_tools.jl")
    include("tle_tools.jl")
    include("data_process.jl")
    include("doppler_tools.jl")
    include("processing_tools.jl")
    include("demos.jl")
    include("benchmarks.jl")

    export Rₑ
    export L1_freq
    export h_parms_tcxo, h_parms_ocxo, h_parms_cesium
    export l1ca_code_length
    export l1ca_db_chipping_rate
    export l1ca_chipping_rate
    export l1ca_codes, define_l1ca_code_type
    export L5_freq
    export L5_code_length
    export nh10_code_length
    export nh20_code_length
    export L5_chipping_rate
    export nh_chipping_rate
    export L5_db_chipping_rate
    export l5q_codes
    export nh10
    export nh20
    export l5i_codes
    export definesignal, definesignal!, definereplica, definereplica!
    export generatesignal!, generatereplica!
    export fft_correlate
    export GNSSData
    export loaddata, reloaddata!
    export gencorrresult
    export courseacquisition, courseacquisition!
    export v_t2p_fa, p_fa2v_t, false_alarm_pdf, detection_pdf, v_t2p_d, p_d2v_t
    export fineacquisition
    export gnsstypes
    export calcinitcodephase
    export calccodeidx
    export calctvector
    export calcsnr
    export loadrinex
    export trackprn
    export calcA, calcC, calcQ
    export plotresults
    export dkalman
    export demo
    export courseacquisitiontest
    export meshgrid
    export data_info_from_name
    export initorbitinfo
    export calcdoppler, calcdopplerrate
    export getGPSSatnums
    export getTLEs
    export GPSData
    export dataprocess
    export doppler2chips, get_chips_and_ϕ
    export define_constellation
    export calcelevation
    export course_acq_est
    export plot_satellite_orbit
    export doppler_distribution
    export definecodetype, definesignaltype
    export get_eop
    export generate_phase_noise, generate_phase_noise!, plot_spectrum
    export calc_code_val
    export process
    export CN0_monte_carlo
    export define_l5_code_type, define_l1ca_code_type
    export get_distribution_bounds
    export binary2OneOrNegativeOne
    export make3dplot

end # module
