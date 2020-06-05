__precompile__()

module GNSSTools

    using Reexport

    using BenchmarkTools

    println("Importing Random...")
    using Random
    println("Importing FFTW...")
    using FFTW
    println("Importing Threads...")
    using Base.Threads
    println("Importing ProgressMeter...")
    using ProgressMeter
    println("Importing Statistics...")
    using Statistics
    println("Importing LinearAlgebra...")
    using LinearAlgebra
    println("Importing SatelliteToolbox...")
    using SatelliteToolbox
    println("Importing PyPlot...")
    using PyPlot
    pygui(true)

    println("Importing module files...")

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

    println("Exporting vars...")

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
    export calcdoppler
    export trackprn
    export calcA, calcC, calcQ
    export plotresults
    export dkalman
    export demo

end # module
