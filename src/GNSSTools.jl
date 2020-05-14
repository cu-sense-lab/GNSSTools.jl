module GNSSTools

    using Reexport

    @reexport using Random
    @reexport using FFTW
    @reexport using Base.Threads
    @reexport using ProgressMeter
    @reexport using Statistics
    @reexport using LinearAlgebra
    @reexport using ControlSystems
    @reexport using SatelliteToolbox
    @reexport using PyPlot
    pygui(true)

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
    include("testGNSSTools.jl")

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
    export demo

end # module
