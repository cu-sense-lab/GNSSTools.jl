"""
    FineAcquisitionResults

A struct that stores the fine acquisition fine acquisition
results for both the carrier and FFT based methods.
"""
struct FineAcquisitionResults{T1}
    prn::Int64
    type::String
    fd_course::Float64
    fd_rate::Float64
    n0_idx_course::Int64
    t_length::Float64
    fd_fine::Float64
    fd_est::Float64
    phi_init::Float64
    M::T1
    P::Array{Float64,2}
    R::Array{Float64,1}
end


"""
    fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                    n₀_idx_course, type::Val{:fft}; fd_rate=0.,
                    t_length=replica.t_length, freq_lim=10000.,
                    σω=10., err_bin_num_ϕ=1, err_bin_num_f=2)

Performs an FFT based fine acquisition on `data`. Note that `t_length` must
equal `replica.t_length`. `data` can be either a `GNSSData` or `L5QSignal`
struct, however, `data` and `replica` must be two seperate structs.
"""
function fineacquisition(data::GNSSSignal, replica::ReplicaSignals, prn, fd_course,
                         n₀_idx_course, type::Val{:fft}; fd_rate=0.,
                         t_length=replica.t_length, freq_lim=10000.,
                         σω=1000., err_bin_num_ϕ=1, err_bin_num_f=2)
    # Generate replica
    # Set signal parameters
    definesignal!(replica;
                  prn=prn, f_d=fd_course,
                  fd_rate=fd_rate, ϕ=0., f_if=0.,
                  # include_carrier=true,
                  # include_thermal_noise=false,
                  # include_phase_noise=false,
                  # include_adc=false,
                  code_start_idx=n₀_idx_course,
                  isreplica=true,
                  noexp=true)
    # Generate signal
    generatesignal!(replica, replica.isreplica)
    # Wipeoff IF and course Doppler from data and multiply by replica
    # sig = data.data.*exp.(-2π.*(data.f_if+fd_course).*data.t.*1im).*replica.data
    @threads for i in 1:replica.sample_num
        # @inbounds replica.data[i] = data.data[i]*exp(-2π*(data.f_if+fd_course)*data.t[i]*1im)*replica.data[i]
        @inbounds replica.data[i] = data.data[i]*cis(-2π*(data.f_if+fd_course+0.5*fd_rate*data.t[i])*data.t[i])*replica.data[i]
        # @inbounds replica.data[i] = data.data[i]*replica.data[i]
    end
    # Perform in place FFT of `replica.data`
    fft!(replica.data)
    # replica.data = fft(sig)
    # Find peak within ±[x]kHz, where `x` is defined by freq_lim
    # From index 1 to N/2: positive frequencies
    # From index N/2 to N: negative frequencies
    Δf = 1/replica.t_length
    N = replica.sample_num
    N_over_2 = Int64(floor(N/2))
    N_lim = Int64(floor(freq_lim/Δf))
    pos_freq_interval = (1, N_lim)
    neg_freq_interval = (N-N_lim, N)
    pk_val = 0. + 0im
    pk_valabs2 = 0.
    pk_idx = 1
    #### TODO: Make the peak searching loops below functions and multi-thread them. ####
    # Check positive frequencies for peak
    @inbounds for i in pos_freq_interval[1]:pos_freq_interval[2]
        data_abs2 = abs2(replica.data[i])
        if data_abs2 > pk_valabs2
            # global pk_valabs2
            # global pk_val
            # global pk_idx
            pk_valabs2 = data_abs2
            pk_val = replica.data[i]
            pk_idx = i
        end
    end
    # Check negative frequencies for peak
    @inbounds for i in neg_freq_interval[1]:neg_freq_interval[2]
        data_abs2 = abs2(replica.data[i])
        if data_abs2 > pk_valabs2
            # global pk_valabs2
            # global pk_val
            # global pk_idx
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
    # Calculate the covariance matrix
    # We estimate the error to be ±2 Doppler bin for the frequency error
    # and ±1 for the phase error.
    pk_low_idx = ((pk_idx-err_bin_num_ϕ)+replica.sample_num)%replica.sample_num
    pk_high_idx = ((pk_idx+err_bin_num_ϕ)+replica.sample_num)%replica.sample_num
    if pk_low_idx == 0
        pk_low_idx = replica.sample_num
    end
    if pk_high_idx == 0
        pk_high_idx = replica.sample_num
    end
    pk_low = replica.data[pk_low_idx]
    pk_high = replica.data[pk_high_idx]
    ϕ_low = atan(imag(pk_low)/real(pk_low))
    ϕ_high = atan(imag(pk_high)/real(pk_high))
    ϕ_init_err = mean([abs(ϕ_init-ϕ_low), abs(ϕ_init-ϕ_high)])
    replica.isreplica = false
    P = diagm([ϕ_init_err^2, (2π*err_bin_num_f*Δf)^2, σω^2])
    R = [ϕ_init_err^2]
    # Return `FineAcquisitionResults` struct
    return FineAcquisitionResults(prn, String(:fft), fd_course, fd_rate, n₀_idx_course,
                                  t_length, fd_fine, fd_est, ϕ_init, "N/A", P, R)
end


"""
    fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                    n₀_idx_course, type::Val{:fft}; fd_rate=0.,
                    t_length=replica.t_length, freq_lim=10000.,
                    σω=10., err_bin_num_ϕ=1, err_bin_num_f=2)

Performs an FFT based fine acquisition on `data`. Note that `t_length` must
equal `replica.t_length`. `data` can be either a `GNSSData` or `L5QSignal`
struct, however, `data` and `replica` must be two seperate structs.
"""
function fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                         n₀_idx_course, type::Val{:fft}; fd_rate=0.,
                         t_length=replica.t_length, freq_lim=10000.,
                         σω=10., err_bin_num_ϕ=1, err_bin_num_f=2)
    # Generate replica
    # Set signal parameters
    definesignal!(replica;
                  prn=prn, f_d=fd_course,
                  fd_rate=fd_rate, ϕ=0., f_if=0.,
                  # include_carrier=true,
                  # include_noise=false,
                  # include_adc=false,
                  code_start_idx=n₀_idx_course,
                  isreplica=true,
                  noexp=true)
    # Generate signal
    generatesignal!(replica, replica.isreplica)
    # Wipeoff IF and course Doppler from data and multiply by replica
    # sig = data.data.*exp.(-2π.*(data.f_if+fd_course).*data.t.*1im).*replica.data
    @threads for i in 1:replica.sample_num
        # @inbounds replica.data[i] = data.data[i]*exp(-2π*(data.f_if+fd_course)*data.t[i]*1im)*replica.data[i]
        @inbounds replica.data[i] = data.data[i]*cis(-2π*(data.f_if+fd_course+0.5*fd_rate*data.t[i])*data.t[i])*replica.data[i]
    end
    # Perform in place FFT of `replica.data`
    fft!(replica.data)
    # replica.data = fft(sig)
    # Find peak within ±[x]kHz, where `x` is defined by freq_lim
    # From index 1 to N/2: positive frequencies
    # From index N/2 to N: negative frequencies
    Δf = 1/replica.t_length
    N = replica.sample_num
    N_over_2 = Int64(floor(N/2))
    N_lim = Int64(floor(freq_lim/Δf))
    pos_freq_interval = (1, N_lim)
    neg_freq_interval = (N-N_lim, N)
    pk_val = 0. + 0im
    pk_valabs2 = 0.
    pk_idx = 1
    #### TODO: Make the peak searching loops below functions and multi-thread them. ####
    # Check positive frequencies for peak
    @inbounds for i in pos_freq_interval[1]:pos_freq_interval[2]
        data_abs2 = abs2(replica.data[i])
        if data_abs2 > pk_valabs2
            # global pk_valabs2
            # global pk_val
            # global pk_idx
            pk_valabs2 = data_abs2
            pk_val = replica.data[i]
            pk_idx = i
        end
    end
    # Check negative frequencies for peak
    @inbounds for i in neg_freq_interval[1]:neg_freq_interval[2]
        data_abs2 = abs2(replica.data[i])
        if data_abs2 > pk_valabs2
            # global pk_valabs2
            # global pk_val
            # global pk_idx
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
    # Calculate the covariance matrix
    # We estimate the error to be ±2 Doppler bin for the frequency error
    # and ±1 for the phase error.
    pk_low_idx = ((pk_idx-err_bin_num_ϕ)+replica.sample_num)%replica.sample_num
    pk_high_idx = ((pk_idx+err_bin_num_ϕ)+replica.sample_num)%replica.sample_num
    if pk_low_idx == 0
        pk_low_idx = replica.sample_num
    end
    if pk_high_idx == 0
        pk_high_idx = replica.sample_num
    end
    pk_low = replica.data[pk_low_idx]
    pk_high = replica.data[pk_high_idx]
    ϕ_low = atan(imag(pk_low)/real(pk_low))
    ϕ_high = atan(imag(pk_high)/real(pk_high))
    ϕ_init_err = mean([abs(ϕ_init-ϕ_low), abs(ϕ_init-ϕ_high)])
    replica.isreplica = false
    P = diagm([ϕ_init_err^2, (2π*err_bin_num_f*Δf)^2, σω^2])
    R = [ϕ_init_err^2]
    # Return `FineAcquisitionResults` struct
    return FineAcquisitionResults(prn, String(:fft), fd_course, fd_rate, n₀_idx_course,
                                  t_length, fd_fine, fd_est, ϕ_init, "N/A", P, R)
end


"""
    fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                    n₀_idx_course, type::Val{:carrier}; fd_rate=0.,
                    t_length=replica.t_length, freq_lim=50000., M=1)

Performs an carrier based fine acquisition on `data`.
`replica` decides what signal type to use and the length of each `M` segment.

`M` is the multiple of the length of replica to be used for fine acquisition.
"""
function fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                         n₀_idx_course, type::Val{:carrier}; fd_rate=0.,
                         t_length=replica.t_length, freq_lim=50000., M=1)
    # Check that the data length is at least 2x greater than the size of `replica`
    Mblocks = Int64(floor(data.sample_num)/(M*replica.sample_num))
    if Mblocks < 2
        errmsg = string("Cannot perform carrier based fine acquisition on data",
                        "\n`data.t_length` must be at least ≳ 2(`replica.t_length`).")
        error(errmsg)
    end
    # Samples per `M` data segment
    N = replica.sample_num
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
        ts = view(data.t, (i-1)*N+1:i*N)
        # Wipeoff IF and course Doppler from data
        @threads for i in 1:replica.sample_num
            @inbounds data.data[i]*exp(-2π*(data.f_if+fd_course)*data.t[i]*1im)*replica.data[i]
        end
        # Perform in place fft operation
        fft!(replica.data)
        # Get FFT peak value
        pk = maximum(abs2.(replica.data))
        # Take square root of `pk`
        pk = sqrt(pk)
        # Calculate ϕ for current processing block
        @inbounds ϕ[i] = atan(imag(pk), real(pk))
        if i == 1
            @inbounds ϕ_init = atan(imag(pk)/real(pk))
        end
    end
    # Take the difference between ϕs & check and correct phase values if > than π
    dϕ = diff(ϕ)
    for i in 1:Mblocks-1
        if dϕ[i] > π
            dϕ[i] -= 2π
        elseif dϕ[i] < π
            dϕ[i] += 2π
        else
            # Do nothing. Within bounds.
        end
    end
    # Compute variance from average phase change
    dϕvar = var(dϕ)
    # Find the maxumum variance and omit from fine Doppler fequency calculation
    maxvaridx = argmax(dϕvar)
    # Take the mean while ommitting the largest phase change in `dϕ`
    dϕavg = 0.
    for i in 1:Mblocks-1
        if i != maxvaridx
            dϕavg += dϕ[i]
        end
    end
    dϕavg /= (Mblocks-2)
    # Calculate the fine Doppler frequency
    fd_fine = dϕavg/(2π*N/f_s)
    fd_est = fd_course + fd_fine
    definesignal!(replica; isreplica=false)
    # Return `FineAcquisitionResults` struct
    return FineAcquisitionResults(prn, String(:carrier), fd_course, fd_rate, n₀_idx_course,
                                  t_length, fd_fine, fd_est, ϕ_init, M)
end
