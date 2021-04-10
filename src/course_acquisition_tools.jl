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
                      prn; fd_center=0., fd_range=5000.,
                      fd_rate=0., Δfd=1/replica.t_length,
                      threads=nthreads(), operation="replace",
                      start_idx=1)


Initializes `corr_result` array and runs `courseacquisition!` method. Returns
`corr_result` and the course acquired code phase, Doppler frequency and SNR.


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
                           prn; fd_center=0., fd_range=5000.,
                           fd_rate=0., Δfd=1/replica.t_length,
                           start_idx=1, return_corrresult=false)
    # Allocate space for correlation result
    corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
    # Perform course acquisition
    courseacquisition!(corr_result, data, replica, prn;
                       fd_center=fd_center, fd_range=fd_range,
                       fd_rate=fd_rate, Δfd=Δfd, start_idx=start_idx)
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
                  exp.(-2π.*data.f_if.*data.t[1:dsize].*1im))
    # Number of bits representing `data`
    nADC = data.nADC
    # Number of Doppler bins
    doppler_bin_num = floor(Int, 2*fd_range/Δfd+1)
    @inbounds for i in 1:doppler_bin_num
        # Calculate Doppler frequency for `i` Doppler bin
        f_d = (fd_center-fd_range) + (i-1)*Δfd
        # Set signal parameters
        definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, phi=0., f_if=0.,
                      code_start_idx=1,
                      nADC=nADC, isreplica=true,
                      noexp=false)
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
    # noise_val = sum(corr_result[max_idx[1],:])
    # PN = (noise_val - pk_val)/(M - 1)
    noise_val = sum(corr_result)
    PN = (noise_val - pk_val)/(N*M - 1)
    SNR_est = 10*log10(sqrt(PS/PN))
    return (n0_est, fd_est, SNR_est)
end
