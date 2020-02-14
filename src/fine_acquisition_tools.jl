"""
    FineAcquisitionResults

A struct that stores the fine acquisition fine acquisition
results for both the carrier and FFT based methods.
"""
struct FineAcquisitionResults{T}
    type::String
    fd_course::Float64
    fd_rate::Float64
    n0_idx_course::Float64
    t_length::Float64
    fd_fine::Float64
    fd_est::Float64
    ϕ_init::Float64
    M::T
end


"""
    fineacquisition(data::GNSSData, replica, fd_course,
                    n₀_course, type::Val{:fft})

Performs an FFT based fine acquisition on `data`. Note that `t_length` must
equal `replica.t_length`.
"""
function fineacquisition(data::GNSSData, replica, fd_course,
                         n₀_idx_course, type::Val{:fft}; fd_rate=0.,
                         t_length=replica.t_length, freq_lim=50000.)
    # Generate replica
    # Set signal parameters
    definesignal!(replica;
                  prn=prn, f_d=fd_course,
                  fd_rate=fd_rate, ϕ=0., f_if=0.,
                  include_carrier=true,
                  include_noise=false,
                  include_adc=false,
                  code_start_idx=n₀_idx_course,
                  isreplica=true)
    # Generate signal
    generatesignal!(replica)
    # Wipeoff IF and course Doppler from data and multiply by replica
    @threads for i in 1:replica.sample_num
        @inbounds data.data[i]*exp(-2π*(data.f_if+fd_course)*data.t[i]*1im)*replica.data[i]
    end
    # Perform in place FFT of `replica.data`
    fft!(replica.data)
    # Find peak within ±[x]kHz, where `x` is defined by freq_lim
    # From index 1 to N/2: positive frequencies
    # From index N/2 to N: negative frequencies
    Δf = 1/data.t_length
    N = data.sample_num
    N_over_2 = Int64(floor(N/2))
    N_lim = Int64(floor(freq_lim/Δf))
    pos_freq_interval = (1, N_lim)
    neg_freq_interval = (N-N_lim, N)
    pk_val = 0. + 0im
    pk_valabs2 = max_val
    pk_idx = 1
    #### TODO: Make the peak searching loops below functions and try to multi-thread them. ####
    # Check positive frequencies for peak
    @inbounds for i in pos_freq_interval[1]:pos_freq_interval[2]
        data_abs2 = abs2(replica.data[i])
        if data_abs2 > pk_valabs2
            pk_valabs2 = data_abs2
            pk_val = replica.data[i]
            pk_idx = i
        end
    end
    # Check negative frequencies for peak
    @inbounds for i in neg_freq_interval[1]:neg_freq_interval[2]
        data_abs2 = abs2(replica.data[i])
        if data_abs2 > pk_valabs2
            pk_valabs2 = data_abs2
            pk_val = replica.data[i]
            pk_idx = i
        end
    end
    # Calculate peak frequency (the Doppler frequency error)
    if pk_idx > N_over_2
        fd_fine = (pk_idx-N-1)*Δf
    else
        fd_fine = (pk_idx-1)*Δf
    end
    fd_est = fd_course + fd_fine
    # Calculate initial phase
    ϕ_init = atan(imag(pk_val)/real(pk_val))
    # Return `FineAcquisitionResults` struct
    return FineAcquisitionResults(String(type), fd_course, fd_rate, n0_idx_course,
                                  t_length, fd_fine, fd_est, ϕ_init)
end


"""
    fineacquisition(data::GNSSData, replica, fd_course,
                    n₀_course, type::Val{:carrier})

Performs an carrier based fine acquisition on `data`.
`replica` decides what signal type to use and the length of each `M` segment.

`M` is the multiple of the length of replica to be used for fine acquisition.
"""
function fineacquisition(data::GNSSData, replica, fd_course,
                         n₀_idx_course, type::Val{:carrier}; fd_rate=0.,
                         t_length=replica.t_length, freq_lim=50000., M=1)
    # Check that the data length is at least 2x greater than the size of `replica`
    Mblocks = Int64(floor(data.sample_num)/(M*replica.sample_num))
    if Mblocks < 2
        errmsg = string("ERROR: Cannot perform carrier based fine acquisition on data",
                        "\n`data.t_length` must be at least ≳ 2(`replica.t_length`).")
        error(errmsg)
    end
    # Samples per `M` data segment
    N = replica.smple_num
    # Set ϕ_init and ϕ
    ϕ_init = 0.
    ϕ = Array{Float64}(undef, Mblocks)
    for i in 1:Mblocks
        # Generate replica
        # Set signal parameters
        definesignal!(replica;
                      prn=prn, f_d=fd_course,
                      fd_rate=fd_rate, ϕ=0., f_if=0.,
                      include_carrier=true,
                      include_noise=false,
                      include_adc=false,
                      code_start_idx=n₀_idx_course,
                      isreplica=true)
        # Generate signal
        generatesignal!(replica)
        # Get a `view` of the current data segment and its corresponding time array
        datasegment = view(data.data, (i-1)*N+1:i*N)
        ts = view(data.t, (i-1)*N+1:i*N))
        # Wipeoff IF and course Doppler from data
        @threads for i in 1:replica.sample_num
            @inbounds data.data[i]*exp(-2π*(data.f_if+fd_course)*data.t[i]*1im)*replica.data[i]
        end
        # Perform in place fft operation
        fft!(replca.data)
        # Get FFT peak value
        pk = maximum(replica.data)
        # Calculate ϕ for current processing block
        @inbounds ϕ[i] = atan(imag(pk), real(pk))
        if i == 1
            @inbounds ϕ_init = atan(imag(pk)/real(pk))
        end
    end
    # 
end
