module GNSSTools

    using Reexport

    @reexport using Random
    @reexport using FFTW
    @reexport using PyPlot
    @reexport using Base.Threads
    @reexport using ProgressMeter
    pygui(true)

    include("constants.jl")
    include("signal_generator.jl")
    include("l5q_code_generator.jl")
    include("gnss_tools.jl")

    export c
    export k
    export L5_freq
    export L5_code_length
    export nh_code_length
    export L5_chipping_rate
    export nh_chipping_rate
    export L5QSignal
    export definesignal, definesignal!
    export generatesignal!
    export fft_correlate
    export GNSSData
    export loaddata
    export gencorrresult
    export courseacquisition!
    export calcsnr
    export testcourseacquisition
    export testcourseacquisitiondata
    export testnoncoherentintegration

end # module
