"""
    FineAcquisitionResults


A struct that stores the fine acquisition fine acquisition
results for both the carrier and FFT based methods.


Fields:

- `prn::Int64`: PRN number processed
- `type::String`: fine acquisition method used (either `"fft"` or `"carrier"`)
- `fd_course::Float64`: course acquired Doppler frequency in Hz
- `fd_rate::Float64`: signal Doppler rate in Hz/s
- `n0_idx_course::Int64`: course acquired code start index
- `t_length::Float64`: coherent integration time in seconds
- `fd_fine::Float64`: fine Doppler frequency in Hz
- `fd_est::Float64`: (fine + course) Doppler frequency in Hz
- `phi_init::Float64`: initial phase estimate in rads
- `M::T1`: number of integrations performed (carrier method only)
- `P::Array{Float64,2}`: 3x3 diagonal matrix containin initial uncertainties
                         of the phase, Doppler, and Doppler rate estimates
- `R::Array{Float64,1}`: measurement uncertainty in rads
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
                    t_length=replica.t_length,
                    σω=10., err_bin_num_ϕ=1, err_bin_num_f=2,
                    M=1)


Performs an FFT based fine acquisition on `data`. Note that `t_length` must
equal `replica.t_length`. `data` can be either a `GNSSData` or `L5QSignal`
struct, however, `data` and `replica` must be two seperate structs.


Required Arguments:

- `data::GNSSSignal`: either `GNSSData` or `ReplicaSignal` struct
- `replica::ReplicaSignal`: struct to use for replica signal generation
    * `replica.t_length` determines the integration time, `T`, used in FFT fine
      acquisition
- `prn::Int`: PRN number to track
- `fd_course`: course acquired Doppler frequency in Hz
- `n₀_idx_course`: initial code start index location
- `type::Val{:fft}`: `Val(:fft)` uses the FFT based fine acquisition method


Optional Arguments:

- `fd_rate`: the Doppler frequency rate in Hz/s `(default = 0Hz)`
- `σω`: expected uncertainty of Doppler rate in Hz/s `(default = 1000Hz/s)`
- `err_bin_num_ϕ`: number of bins to to look around peak to estimate the phase
                   uncertainty `(default = 1)`
- `err_bin_num_f`: number of frequency bins equal to the frequency uncertainty
                   `(default = 2)`
    * filter σ is equal to `err_bin_num_f×(1/t_length)`
- `M`: the number of consecutive integrations and use the one with the highest SNR
    * `(default = 1)`


Returns:

- `FineAcquisitionResults` struct
"""
function fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                         n₀_idx_course, type::Val{:fft}; fd_rate=0.,
                         σω=1000., err_bin_num_ϕ=1, err_bin_num_f=0.25,
                         M=1)
    t_length=replica.t_length
    Δf = 1/replica.t_length
    N = replica.sample_num
    N_over_2 = Int64(floor(N/2))
    pk_vals = Array{Complex{Float64}}(undef, M)
    pk_valabs2s = Array{Float64}(undef, M)
    pk_idxs = Array{Int}(undef, M) 
    PNs = Array{Float64}(undef, M)
    SNR_ests = Array{Float64}(undef, M)
    fd_fines = Array{Float64}(undef, M)
    fd_ests = Array{Float64}(undef, M)
    ϕ_inits = Array{Float64}(undef, M)
    ϕ_init_errs = Array{Float64}(undef, M)
    for m in 1:M
        # Generate replica
        # Set signal parameters
        definereplica!(replica;
                    prn=prn, f_d=fd_course,
                    fd_rate=fd_rate, phi=0., f_if=0.,
                    code_start_idx=n₀_idx_course,
                    include_carrier=false,
                    start_t=(m-1)*t_length)
        # Generate signal
        generatereplica!(replica)
        # Wipeoff IF and course Doppler from data and multiply by replica
        @threads for i in 1:replica.sample_num
            j = (m-1)*N + i
            t = calc_t_at_i(j, (m-1)*t_length, data.f_s)
            @inbounds replica.data[i] = data.data[j]*cis(-2π*(data.f_if+fd_course+0.5*fd_rate*t)*t)*replica.data[i]
        end
        # Perform in place FFT of `replica.data`
        fft!(replica.data)
        # Find peak within ±[x]kHz, where `x` is defined by freq_lim
        # From index 1 to N/2: positive frequencies
        # From index N/2 to N: negative frequencies
        # N_lim = Int64(floor(freq_lim/Δf))
        # pos_freq_interval = (1, N_lim)
        # neg_freq_interval = (N-N_lim, N)
        pk_val = 0. + 0im
        pk_valabs2 = 0.
        pk_idx = 1
        PN = 0.
        #### TODO: Make the peak searching loops below functions and multi-thread them. ####
        # Check positive frequencies for peak
        @inbounds for i in 1:N
            data_abs2 = abs2(replica.data[i])
            PN += data_abs2
            if data_abs2 > pk_valabs2
                pk_valabs2 = data_abs2
                pk_val = replica.data[i]
                pk_idx = i
            end
        end
        PN = (PN - pk_valabs2) / (N - 1)
        SNR_est = calc_snr(pk_valabs2, PN)
        # # Check negative frequencies for peak
        # @inbounds for i in neg_freq_interval[1]:neg_freq_interval[2]
        #     data_abs2 = abs2(replica.data[i])
        #     if data_abs2 > pk_valabs2
        #         pk_valabs2 = data_abs2
        #         pk_val = replica.data[i]
        #         pk_idx = i
        #     end
        # end
        # Calculate peak frequency (the Doppler frequency error)
        if pk_idx > N_over_2
            fd_fine = (pk_idx-N-1)*Δf
        else
            fd_fine = (pk_idx-1)*Δf
        end
        fd_est = fd_course + fd_fine
        # Calculate initial phase
        ϕ_init = atan(imag(pk_val), real(pk_val))
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

        pk_vals[m] = pk_val
        pk_valabs2s[m] = pk_valabs2
        pk_idxs[m] = pk_idx
        PNs[m] = PN
        SNR_ests[m] = SNR_est
        fd_fines[m] = fd_fine
        fd_ests[m] = fd_est
        ϕ_inits[m] = ϕ_init
        ϕ_init_errs[m] = ϕ_init_err 
    end
    idx = argmax(SNR_ests)
    pk_val = pk_vals[idx]
    pk_valabs2 = pk_valabs2s[idx]
    pk_idx = pk_idxs[idx]
    PN = PNs[idx]
    SNR_est = SNR_ests[idx]
    fd_fine = fd_fines[idx]
    fd_est = fd_ests[idx]
    ϕ_init = ϕ_inits[idx]
    ϕ_init_err = ϕ_init_errs[idx]
    # replica.isreplica = false
    P = diagm([ϕ_init_err^2, (2π*err_bin_num_f*Δf)^2, (2π*σω)^2])
    R = [ϕ_init_err^2]
    # Return `FineAcquisitionResults` struct
    return FineAcquisitionResults(prn, String(:fft), fd_course, fd_rate, n₀_idx_course,
                                  t_length, fd_fine, fd_est, ϕ_init, "N/A", P, R)
end


"""
    fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                    n₀_idx_course, type::Val{:carrier}; fd_rate=0.,
                    t_length=replica.t_length, freq_lim=50000., M=10,
                    σω=1000.)


Performs an carrier based fine acquisition on `data`.
`replica` decides what signal type to use and the length of each `M` segment.

`M` is the multiple of the length of replica to be used for fine acquisition.


Required Arguments:

- `data::GNSSSignal`: either `GNSSData` or `ReplicaSignal` struct
- `replica::ReplicaSignal`: struct to use for replica signal generation
    * `replica.t_length` determines the integration time, `T`
- `prn::Int`: PRN number to track
- `fd_course`: course acquired Doppler frequency in Hz
- `n₀_idx_course`: initial code start index location
- `type::Val{:fft}`: `Val(:fft)` uses the FFT based fine acquisition method


Optional Arguments:

- `fd_rate`: the Doppler frequency rate in Hz/s `(default = 0Hz)`
- `t_length`: `[NOT USED]` `(default = replica.t_length)`
- `M`: Number of integrations to perform `(default = 10)`
- `σω`: expected uncertainty of Doppler rate in Hz/s `(default = 1000Hz/s)`


Returns:

- `FineAcquisitionResults` struct
"""
function fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                         n₀_idx_course, type::Val{:carrier}; fd_rate=0.,
                         t_length=replica.t_length, M=10, σω=1000.)
    # Check that there are at least M replica signals that can fit in time
    # into data.
    Mblocks = floor(Int, data.sample_num/replica.sample_num)
    if Mblocks < M
        error("Data must be at least $(M*replica.t_length) seconds long.")
    end
    # Samples per `M` data segment
    N = replica.sample_num
    # Set ϕ_init and ϕ
    ϕ_init = 0.
    ϕ = Array{Float64}(undef, M)
    ϕi = 0.
    for i in 1:M
        # Generate replica
        # Set signal parameters
        definereplica!(replica;
                       prn=prn, f_d=fd_course,
                       fd_rate=fd_rate, phi=0, f_if=0.,
                       code_start_idx=n₀_idx_course,
                       include_carrier=false)
        # Generate signal
        generatereplica!(replica)
        # Get a `view` of the current data segment and its corresponding time array
        datasegment = view(data.data, (i-1)*N+1:i*N)
        # Wipeoff IF and course Doppler from data
        @threads for j in 1:replica.sample_num
            t = calc_t_at_i((i-1)*N + j, 0, data.f_s)
            @inbounds replica.data[j] = datasegment[j]*cis(-2π*(data.f_if+fd_course+0.5*fd_rate*t)*t)*replica.data[j]
        end
        # Perform in place fft operation
        fft!(replica.data)
        # Get FFT peak value
        pk_idx = argmax(abs2.(replica.data))
        # Take square root of `pk`
        pk = replica.data[pk_idx]
        # Calculate ϕ for current processing block
        ϕi = atan(imag(pk), real(pk))
        # ϕi = atan(imag(pk)/real(pk))
        @inbounds ϕ[i] = ϕi
        if i == 1
            @inbounds ϕ_init = ϕi
        end
    end
    # Take the difference between ϕs & check and correct phase values if > than π
    dϕ = diff(ϕ)
    for i in 1:M-1
        if dϕ[i] > π
            dϕ[i] -= 2π
        elseif dϕ[i] < -π
            dϕ[i] += 2π
        else
            # Do nothing. Within bounds.
        end
    end
    # Compute variance from average phase change
    # dϕvar = var(dϕ)
    dϕvar = abs2.(dϕ)
    # Find the maxumum variance and omit from fine Doppler fequency calculation
    maxvaridx = argmax(dϕvar)
    # Take the mean while ommitting the largest phase change in `dϕ`
    dϕavg = 0.
    for i in 1:M-1
        if i != maxvaridx
            dϕavg += dϕ[i]
        end
    end
    dϕavg /= (M-2)
    # Calculate the fine Doppler frequency
    f_s = replica.f_s
    # fd_fine = dϕavg/(2π*N/f_s)
    fd_fine = dϕavg/(2π*replica.t_length)
    fd_est = fd_course + fd_fine
    dϕ_init_err = std(([dϕ[1:maxvaridx-1]; dϕ[maxvaridx+1:end]] .- dϕavg))
    fd_err = dϕ_init_err/(2π*replica.t_length)
    # P = diagm([dϕ_init_err^2, 74.7*fd_err^2, σω^2])
    # P = diagm([dϕ_init_err^2, fd_err^2, σω^2])
    P = diagm([dϕ_init_err^2, (2π*replica.t_length)^2, (2π*σω)^2])
    R = [dϕ_init_err^2]
    # Return `FineAcquisitionResults` struct
    return FineAcquisitionResults(prn, String(:carrier), fd_course, fd_rate, n₀_idx_course,
                                  t_length, fd_fine, fd_est, ϕ_init, "N/A", P, R)
end


"""
    fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                    n₀_idx_course, type::Val{:ffl}; fd_rate=0.,
                    t_length=replica.t_length, freq_lim=50000., M=10,
                    σω=1000.)


Performs an carrier based fine acquisition on `data`.
`replica` decides what signal type to use and the length of each `M` segment.

`M` is the multiple of the length of replica to be used for fine acquisition.


Required Arguments:

- `data::GNSSSignal`: either `GNSSData` or `ReplicaSignal` struct
- `replica::ReplicaSignal`: struct to use for replica signal generation
    * `replica.t_length` determines the integration time, `T`
- `prn::Int`: PRN number to track
- `fd_course`: course acquired Doppler frequency in Hz
- `n₀_idx_course`: initial code start index location
- `type::Val{:fll}`: `Val(:fll)` uses the FLL based fine acquisition method


Optional Arguments:

- `fd_rate`: the Doppler frequency rate in Hz/s `(default = 0Hz)`
- `t_length`: `[NOT USED]` `(default = replica.t_length)`
- `M`: Number of integrations to perform `(default = 10)`
- `σω`: expected uncertainty of Doppler rate in Hz/s `(default = 1000Hz/s)`
- `B`: FLL bandwith in Hz `(default = 20)`
- `dampening`: dampening ratio for FLL `(default = 0.707)`


Returns:

- `FineAcquisitionResults` struct
"""
function fineacquisition(data::GNSSSignal, replica::ReplicaSignal, prn, fd_course,
                         n₀_idx_course, type::Val{:fll}; fd_rate=0.,
                         t_length=replica.t_length, M=10, σω=1000.,
                         B=20, dampening=0.707, q_a=1, h₋₂=2e-20)
    ωₙ = B/0.53
    T = t_length
    sig_freq = replica.signal_type.sig_freq
	qω = 2*π^2*h₋₂  # oscillator frequency PSD
    # Check that there are at least M replica signals that can fit in time
    # into data.
    Mblocks = floor(Int, data.sample_num/replica.sample_num)
    if Mblocks < M
        error("Data must be at least $(M*replica.t_length) seconds long.")
    end
    # Samples per `M` data segment
    N = replica.sample_num
    # Initialize FLL
    Δω = Array{Float64}(undef, M-1)
    SNR = Array{Float64}(undef, M)
    x = Array{Float64}(undef, 2, M)
    x[:,1] .= [fd_course, fd_rate] .* 2π
    P = Array{Float64}(undef, 2, M)
    P[:,1] .= ([0.5/t_length, σω] .* 2π).^2
    x⁺ = [x[1,1]; x[2,1]]
    P⁺ = [P[1,1] 0; 0 P[2,1]]
    Δθₖ₋₁ = 0
    K = [2*dampening*ωₙ*T-3*ωₙ^2*T^2/2; ωₙ^2*T]
    K[2] = 0
    H = [1 T/2]
    A = [1 T; 0 1]
    Q = (2π*sig_freq)^2 .* [T*qω+T^3*q_a/(3*c^2) T^2*q_a/(2*c^2);
                            T^2*q_a/(2*c^2)      T*q_a/c^2]
    for i in 1:M
        # Generate replica
        # Set signal parameters
        f_d, fd_rate = x⁺ ./ 2π
        fd_rate = 0
        definereplica!(replica;
                       prn=prn, f_d=fd_course,
                       fd_rate=fd_rate, phi=0, f_if=0.,
                       code_start_idx=n₀_idx_course,
                       include_carrier=false,
                       start_t=(i-1)*T)
        # Generate signal
        generatereplica!(replica)
        # Get a `view` of the current data segment and its corresponding time array
        datasegment = view(data.data, (i-1)*N+1:i*N)
        # Wipeoff IF and course Doppler from data
        @threads for j in 1:replica.sample_num
            t = calc_t_at_i((i-1)*N + j, 0, data.f_s)
            @inbounds replica.data[j] = datasegment[j]*cis(-2π*(data.f_if+fd_course+0.5*fd_rate*t)*t)*replica.data[j]
        end
        # Perform in place fft operation
        fft!(replica.data)
        # Get FFT peak value
        pk_idx = argmax(abs2.(replica.data))
        # Take square root of `pk`
        pk = replica.data[pk_idx]
        # Calculate SNR
        PS = abs2(pk)
        PN = (sum(abs2.(replica.data)) - PS) / (N - 1)
        snr = calc_snr(PS, PN)
        SNR[i] = snr
        # Calculate ϕ for current processing block
        # ϕi = atan(imag(pk), real(pk))
        Δϕₖ = atan(imag(pk)/real(pk))
        if i == 1
            Δωₖ = 0
        else
            dϕ = (Δϕₖ - Δθₖ₋₁)
            if dϕ > π
                dϕ -= 2π
            elseif dϕ < -π
                dϕ += 2π
            else
                # Do nothing. Within bounds.
            end
            Δωₖ = dϕ/T
            Δω[i-1] = Δωₖ
        end
        Δθₖ₋₁ = Δϕₖ
        # x⁺ = x⁺ + K*Δωₖ
        x⁺ = x⁺ .+ Δωₖ
        # println(Δωₖ/2π)
        P⁺ = (I - K*H)*P⁺
        x[:,i] .= x⁺
        P[:,i] .= diag(P⁺)
        # Propagate states
        x⁺ = A*x⁺
        P⁺ = A*P⁺*A' + Q
    end
    return (x, P, A, Q, H, K, SNR, Δω)
    # Return `FineAcquisitionResults` struct
    # return FineAcquisitionResults(prn, String(:carrier), fd_course, fd_rate, n₀_idx_course,
    #                               t_length, fd_fine, fd_est, ϕ_init, "N/A", P, R)
end
