"""
    gencorrresult(fd_range, Δfd, sample_num; iszeros=false)


Generates a 2D `Array{Float64,2}` to store the cross correlation results output
from the `courseacquisition` function.


Required Arguments:

- `fd_range::Float64`: the range of Doppler plus minus the center frequency in Hz
- `Δfd::Float64`: frequency bin width in Hz


Optional Arguments:

- `iszeros::Bool`: if `true`, sets contents of array to zero, otherwise its
                   elements are set to `undef` `(default = false)`


Returns:

- `Array{Float64,2}`: 2D array to store acquistion result
"""
function gencorrresult(fd_range, Δfd, sample_num; iszeros=false)
    doppler_bin_num = round(Int, 2*fd_range/Δfd+1)
    if iszeros
        # Initialize with zeros
        return zeros(Float64, doppler_bin_num, sample_num)
    else
        # Allocate space, but don't initialize elements to a value
       return Array{Float64}(undef, doppler_bin_num, sample_num)
    end
end


"""
    courseacquisition(data::GNSSSignal, replica::ReplicaSignal,
                      prn; fd_center=0., fd_range=5000., M=1,
                      fd_rate=0., Δfd=1/replica.t_length,
                      threads=nthreads(), operation="replace",
                      start_idx=1)


Performs course acquisition using an integration time specified by
`replica.t_length`, and adds `M` coherent integrations non-coherently.


Required Arguments:

- `data::GNSSSignal`: either a `GNSSData` or `ReplicaSignal`
    * can be either simulated or real data
- `replica::ReplicaSignal`: replica signal whose length determines the
                             integration time used for course acquisition
- `prn::Int`: PRN number refering to code to search for


Optional Arguments:

- `fd_center::Float64`: center Doppler frequency bin value in Hz `(default = 0Hz)`
- `fd_range::Float64`: the range of Doppler frequencies in Hz to search plus and
                       minus the center Doppler frequency bin `(default = 5000Hz)`
- `M::Int`: the number of coherent integrations to add together. This is
            non-coherent integration `(default = 1)`
- `fd_rate::Float64`: the Doppler frequency rate in Hz/s `(default = 0Hz)`
- `Δfd::Float64`: the Doppler bin width in Hz `(default = 1/replica.t_length)`
    * usually determined by the length of the integration time
    * integration time is determined by the length of the replica
      `replica.t_length`
- `start_idx::Int`: starting sample in data plus `replica.sample_num` that is
                    used for processing `(default = 1)`
- `return_corrresult::Bool`: set to `true` to return the 2D correlation result
                             `(default = true)`


Returns:

- `fd_est::Float64`: Doppler frequency bin corresponding to peak in Hz
- `n0_est::Int`: index in data where code starts
- `SNR_est::Float64`: SNR of the correlation peak in dB
- `corr_result::Array{Float64,2}`: (optional) 2D array to store acquistion result
"""
function courseacquisition(data::GNSSSignal, replica::ReplicaSignal,
                           prn; fd_center=0., fd_range=5000., M=1,
                           fd_rate=0., Δfd=1/replica.t_length,
                           start_idx=1, return_corrresult=false)
    # Perform coherent intefration where the integration time is decided based
    # off the length of `replica.t_length`
    if M == 1
        # Allocate space for correlation result
        corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
        # Perform course acquisition
        courseacquisition!(corr_result, data, replica, prn;
                           fd_center=fd_center, fd_range=fd_range,
                           fd_rate=fd_rate, Δfd=Δfd, start_idx=start_idx,
                           operation="replace")
    # Perform non-coherent integration where the integration time is based off
    # `replica.t_length` and `M` coherent integrations are added together
    else
        @assert M > 1
        # Allocate space for correlation result
        corr_result = gencorrresult(fd_range, Δfd, replica.sample_num, 
                                    iszeros=true)
        # Perform course acquisition
        for i in 1:M
            idx = start_idx + (i-1)*replica.sample_num
            replica.t .+= (i > 1)*replica.t_length
            courseacquisition!(corr_result, data, replica, prn;
                               fd_center=fd_center, fd_range=fd_range,
                               fd_rate=fd_rate, Δfd=Δfd, start_idx=idx,
                               operation="add")
        end
        replica.t .= calctvector(replica.sample_num, replica.f_s)
    end
    # Get the code offset in samples, course Doppler, and peak SNR in dB
    n0_est, fd_est, SNR_est = course_acq_est(corr_result, fd_center, fd_range,
                                             Δfd)
    if return_corrresult
        return (fd_est, n0_est, SNR_est, corr_result)
    else
        return (fd_est, n0_est, SNR_est)
    end
end


"""
    courseacquisition!(corr_result::Array{Float64,2},
                       data::GNSSSignal, replica::ReplicaSignal,
                       prn; fd_center=0., fd_range=5000.,
                       fd_rate=0., Δfd=1/replica.t_length,
                       operation="replace", start_idx=1)


Performs course acquisition on either `GNSSData` or `ReplicaSignal`
type struct using defined `ReplicaSignal` type struct. No need
to use `generatereplica!` before calling this function.
Operates in place on `corr_result`.

*NOTE:* If `data` is a replica signal structure type,
        it's structure type must be the same as the type
        of `replica`. They must also be separately defined
        or one must be a deep copy of the other. Do not
        pass the same structure to the `data` and `replica`
        arguments.

`corr_result` contains |conj(fft(replica)*fft(data)|² per
Doppler bin.


Required Arguments:

- `corr_result::Array{Float64,2}`: 2D array to store acquistion result
- `data::GNSSSignal`: either a `GNSSData` or `ReplicaSignal`
    * can be either simulated or real data
- `replica::ReplicaSignal`: replica signal whose length determines the
                             integration time used for course acquisition
- `prn::Int`: PRN number refering to code to search for


Optional Arguments:

- `fd_center::Float64`: center Doppler frequency bin value in Hz `(default = 0Hz)`
- `fd_range::Float64`: the range of Doppler frequencies in Hz to search plus and
                       minus the center Doppler frequency bin `(default = 5000Hz)`
- `fd_rate::Float64`: the Doppler frequency rate in Hz/s `(default = 0Hz)`
- `Δfd::Float64`: the Doppler bin width in Hz `(default = 1/replica.t_length)`
    * usually determined by the length of the integration time
    * integration time is determined by the length of the replica
      `replica.t_length`
- `start_idx::Int`: starting sample in data plus `replica.sample_num` that is
                    used for processing `(default = 1)`



Modifies and Returns:

- `corr_result::Array{Float64,2}`: 2D array to store acquistion result
"""
function courseacquisition!(corr_result::Array{Float64,2},
                            data::GNSSSignal, replica::ReplicaSignal,
                            prn; fd_center=0., fd_range=5000.,
                            fd_rate=0., Δfd=1/replica.t_length,
                            operation="replace", start_idx=1)
    # Number of data samples
    dsize = replica.sample_num
    # Pre-plan FFTs and IFFTs
    pfft = plan_fft!(replica.data)  # In-place FFT plan
    pifft = plan_ifft!(replica.data) # In-place IFFT plan
    # Carrier wipe data signal, make copy, and take FFT
    datafft = fft(data.data[start_idx:start_idx+dsize-1] .*
                  exp.(-2π.*data.f_if.*replica.t[1:dsize].*1im))
    # Number of bits representing `data`
    nADC = data.nADC
    # Number of Doppler bins
    doppler_bin_num = floor(Int, 2*fd_range/Δfd+1)
    @inbounds for i in 1:doppler_bin_num
        # Calculate Doppler frequency for `i` Doppler bin
        f_d = (fd_center-fd_range) + (i-1)*Δfd
        # Set signal parameters
        definereplica!(replica;
                       prn=prn, f_d=f_d,
                       fd_rate=fd_rate, phi=0., f_if=0.,
                       code_start_idx=1,
                       include_carrier=true)
        # Generate signal
        generatereplica!(replica)
        # Perform in place FFT on replica
        pfft*replica.data
        # Take conjugate of FFT(replica) and multiply with FFT of
        # data. The result is stored in `replica.data`
        conjAmultB1D!(replica.data, datafft, dsize)
        # Take IFFT in place
        pifft*replica.data
        # Take `abs2` in place and save result
        # Either replace, add to, or multiply by pre-existing value
        # in `corr_result`
        @threads for j in 1:dsize
            if operation == "replace"
                @inbounds corr_result[i,j] = abs2(replica.data[j])
            elseif operation == "add"
                @inbounds corr_result[i,j] += abs2(replica.data[j])
            elseif operation == "multiply"
                @inbounds corr_result[i,j] *= abs2(replica.data[j])
            end
        end
    end
    # Set `isreplica` flag to false in `replica`
    definesignal!(replica, isreplica=false)
    return corr_result
end


"""
    course_acq_est(corr_result)


Estimate the code start index, Doppler frequency and peak SNR from the
correlation result.


Required Arguments:

- `corr_result::Array{Float64,2}`: 2D array that contains course acquistion
                                   result


Returns:

- `n0_est::Int`: index in data where code starts
- `fd_est::Float64`: Doppler frequency bin corresponding to peak in Hz
- `SNR_est::Float64`: SNR of the peak in dB
"""
function course_acq_est(corr_result, fd_center, fd_range, Δfd)
    # Get peak maximum index location
    max_idx = argmax(corr_result)
    # Calculate course Doppler frequency
    fd_est = (fd_center-fd_range) + (max_idx[1]-1)*Δfd
    # Get course peak index location in time
    n0_est = max_idx[2]
    pk_val = corr_result[max_idx]
    # Compute course acquisition SNR
    PS = pk_val
    N, M = size(corr_result)
    noise_val = sum(corr_result[max_idx[1],:])
    PN = (noise_val - pk_val)/(M - 1)
    # SNR_est = 10*log10(sqrt(PS/PN))
    SNR_est = 10*log10(PS/PN)
    return (n0_est, fd_est, SNR_est)
end


"""
    v_t2p_fa(v_t, σ_n)


Calculate the probability of false alarm based off the threshold and noise
standard deviation.


Required Arguments:

- `v_t`: the threshold set for acquisition
- `σ_n`: the standard deviation of the noise


Returns:

- `p_fa`: the probability of false alarm
"""
function v_t2p_fa(v_t, σ_n)
    return exp(-v_t^2 / (2*σ_n^2))
end


"""
    p_fa2v_t(p_fa, σ_n)


Calculate the threshold value for acquisition based off the probability of 
false alarm and the noise standard deviation.


Required Arguments:

- `p_fa`: the probability of false alarm between 0 and 1
- `σ_n`: the standard deviation of the noise


Returns:

- `v_t`: the threshold for acquisition
"""
function p_fa2v_t(p_fa, σ_n)
    return sqrt(2*σ_n^2*log(1/p_fa))
end


"""
    false_alarm_pdf(v, σ_n)


Calculate the value of the false alarm PDF at a given value, `v`, with
noise standard deviation, `σ_n`.


Required Arguments:

- `v`: value to evaluate PDF at
- `σ_n`: standard deviation of the noise


Returns:

- PDF value at `v`
"""
function false_alarm_pdf(v, σ_n)
    return v*exp(-v^2/σ_n)/σ_n^2
end


"""
    I₀(x)


Evaluates the Modified Bessel function of the 1st kind order zero.


Arguments:

- `x`: where to evaluate I₀ at


Returns:

- I₀(x)
"""
function I₀(x)
    return quadgk(θ -> exp(big(x)*cos(big(θ))), 0, π)[1]/π
end


"""
    detection_pdf(v, σ_n; SNR=missing, v_i=missing)


Calculate the value of the detection PDF at a given value, `v`, with
noise standard deviation, `σ_n`, and `SNR`.


Required Arguments:

- `v`: value to evaluate PDF at
- `σ_n`: standard deviation of the noise
- `SNR`: desired signal-to-noise ratio in dB '(default = missing)` or
- 'v_i': the amplitude of the signature in place of the `SNR` '(default = missing)`



Returns:

- PDF value at `v`
"""
function detection_pdf(v, σ_n; SNR=missing, v_i=missing)
    if ismissing(SNR) && !ismissing(v_i)
        x = v*v_i/σ_n^2
        return (v/σ_n^2)*exp(-(v^2 + v_i^2)/(2*σ_n^2))*I₀(x)
    elseif !ismissing(SNR) && ismissing(v_i)
        SNR_linear = 10^(SNR/10)
        x = (v/σ_n)*sqrt(2*SNR_linear)
        return (v/σ_n^2)*exp(-(v^2/(2*σ_n^2) + SNR_linear))*I₀(x)
    else
        error("Must provide either `SNR` or `v_i`, not both or none.")
    end
end


"""
    v_t2p_d(v_t, σ_n, SNR; rtol=1e-3, order=3, std_scalar=20, 
            steps=missing)


Calculates the probability of detection based off a threshold value, `v_t`,
noise standard deviation, `σ_n`, and desired `SNR`.


Required Arguments:

- `v_t`: threshold value used for acquisition
- `σ_n`: standard deviation of the noise
- `SNR`: desired signal-to-noise ratio in dB '(default = missing)` or
- `v_i`: the amplitude of the signature in place of the `SNR` '(default = missing)`


Optional Arguments:

- `rtol`: tolerance for numerical integration `(default = 1e-3)`
- `order`: order of the numerical integration `(default = 3)`


Returns:

- probability of detection, `P_d`
"""
function v_t2p_d(v_t, σ_n; SNR=missing, v_i=missing, 
                 rtol=1e-10, order=7)
    v_t = big(v_t)
    σ_n = big(σ_n)
    if !ismissing(SNR) && ismissing(v_i)
        SNR = big(SNR)
    elseif ismissing(SNR) && !ismissing(v_i)
        v_i = big(v_i)
    else
        error("Must provide either SNR or v_i.")
    end
    p_d_pdf(v) = detection_pdf(v, σ_n; SNR=SNR, v_i=v_i)
    return quadgk(p_d_pdf, v_t, Inf, rtol=rtol, order=order)[1]
end


"""
    p_d2v_t(p_d, σ_n, SNR)


Calculates the threshold value, `v_t`, based off the probability of detection,
`p_d`, the noise standard deviation, `σ_n`, and the desired `SNR`.


Required Arguments:

- `p_d`: probability of detection between 0 and 1
- `σ_n`: standard deviation of the noise
- `SNR`: desired signal-to-noise ratio in dB


Optional Arguments:

- `rtol`: tolerance for numerical integration `(default = 1e-3)`
- `order`: order of the numerical integration `(default = 3)` 


Returns:

- threshold value, `v_t`
"""
function p_d2v_t(p_d, σ_n, SNR; rtol=1e-10, order=7)
    SNR_linear = 10^(SNR/10)
    p_d_pdf(v) = detection_pdf(v, σ_n; SNR=SNR_linear)
    Pd(v) = quadgk(p_d_pdf, v[1], Inf, rtol=rtol, order=order)[1]
    δp_d(v) = abs(Pd(v) - p_d)
    result = optimize(δp_d, [0.], [1000], [0.])
    return result
end
