__precompile__()

module GNSSTools

    include("greeting.jl")
    export logo, greeting
    printstyled(logo, color=:default, bold=false)
    println(greeting)

    using Reexport

    print("Importing Random...")
    using Random
    print("Done\nImporting FFTW...")
    using FFTW
    print("Done\nImporting Threads...")
    using Base.Threads
    print("Done\nImporting ProgressMeter...")
    using ProgressMeter
    print("Done\nImporting Statistics...")
    using Statistics
    print("Done\nImporting LinearAlgebra...")
    using LinearAlgebra
    print("Done\nImporting CPUTime...")
    using CPUTime
    print("Done\nImporting SatelliteToolbox...")
    using SatelliteToolbox
    print("Done\nImporting ODE...")
    using ODE
    print("Done\nImporting PyPlot...")
    using PyPlot
    Axes3D = PyPlot.Axes3D
    pygui(true)
    println("Done")

    print("Importing module files...")

    include("constants.jl")
    include("general_tools.jl")
    include("signal_types.jl")
    include("calc_code_val.jl")
    include("definesignal.jl")
    include("generatesignal.jl")
    include("l1ca_code_generator.jl")
    include("l5_code_generator.jl")
    include("gnss_data_tools.jl")
    include("course_acquisition_tools.jl")
    include("fine_acquisition_tools.jl")
    include("rinex_tools.jl")
    include("tle_tools.jl")
    include("tracking_tools.jl")
    include("dkalman.jl")
    include("testGNSSTools.jl")

    print("Done\nExporting variables and methods...")

    export c
    export k
    export Rₑ
    export L1_freq
    export l1ca_code_length
    export l1ca_db_chipping_rate
    export l1ca_chipping_rate
    export l1ca_codes
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
    export definesignal
    export definesignal!
    export generatesignal!
    export fft_correlate
    export GNSSData
    export loaddata
    export gencorrresult
    export courseacquisition!
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

    println("Done")

end # module
