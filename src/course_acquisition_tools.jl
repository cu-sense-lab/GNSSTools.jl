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
    if iszeros
        return zeros(Float64, Int(2*fd_range/Δfd+1), sample_num)
    else
       return Array{Float64}(undef, Int(2*fd_range/Δfd+1), sample_num)
    end
end


"""
    courseacquisition(data::GNSSSignal, replica::ReplicaSignals,
                      prn; fd_center=0., fd_range=5000.,
                      fd_rate=0., Δfd=1/replica.t_length,
                      threads=nthreads(), operation="replace",
                      start_idx=1)


Initializes `corr_result` array and runs `courseacquisition!` method. Returns
`corr_result` and the course acquired code phase, Doppler frequency and SNR.


Required Arguments:

- `data::GNSSSignal`:
- `replica::ReplicaSignals`:
- `prn::Int`:


Optional Arguments:

- `fd_center::Float64`:
- `fd_range::Float64`:
- `fd_rate::Float64`:
- `Δfd::Float64`:
- `start_idx::Int`:


Returns:

- `corr_result::Array{Float64,2}`: 2D array to store acquistion result
- `fd_est::Float64`: Doppler frequency bin corresponding to peak in Hz
- `n0_est::Int`: index in data where code starts
- `SNR_est::Float64`: SNR of the correlation peak in dB
"""
function courseacquisition(data::GNSSSignal, replica::ReplicaSignals,
                           prn; fd_center=0., fd_range=5000.,
                           fd_rate=0., Δfd=1/replica.t_length,
                           start_idx=1)
    # Allocate space for correlation result
    corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
    # Perform course acquisition
    courseacquisition!(corr_result, data, replica, prn;
                       fd_center=fd_center, fd_range=fd_range,
                       fd_rate=fd_rate, Δfd=Δfd, start_idx=start_idx)
    n0_est, fd_est, SNR_est = course_acq_est(corr_result)
    return (corr_result, fd_est, n0_est, SNR_est)
end


"""
    courseacquisition!(corr_result::Array{Float64,2},
                       data::GNSSSignal, replica::ReplicaSignals,
                       prn; fd_center=0., fd_range=5000.,
                       fd_rate=0., Δfd=1/replica.t_length,
                       operation="replace", start_idx=1)


Performs course acquisition on either `GNSSData` or `ReplicaSignals`
type struct using defined `ReplicaSignals` type struct. No need
to use `generatesignal!` before calling this function.
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

- `corr_result::Array{Float64,2}`:
- `data::Vector`:
- `replica::ReplicaSignals`:
- `prn::Int`:


Optional Arguments:

- `fd_center::Float64`:
- `fd_range::Float64`:
- `fd_rate::Float64`:
- `Δfd::Float64`:
- `start_idx::Int`:


Modifies and Returns:

- `corr_result::Array{Float64,2}`: 2D array to store acquistion result
"""
function courseacquisition!(corr_result::Array{Float64,2},
                            data::GNSSSignal, replica::ReplicaSignals,
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
    # datafft = fft(data.data)
    # Number of bits representing `data`
    nADC = data.nADC
    # Number of Doppler bins
    doppler_bin_num = Int(fd_range/Δfd*2+1)
    @inbounds for i in 1:doppler_bin_num
        # Calculate Doppler frequency for `i` Doppler bin
        f_d = (fd_center-fd_range) + (i-1)*Δfd
        # Set signal parameters
        definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, phi=0., f_if=0.,
                      # include_carrier=true,
                      # include_thermal_noise=false,
                      # include_phase_noise=false,
                      # include_adc=false,
                      code_start_idx=1,
                      nADC=nADC, isreplica=true,
                      noexp=false)
        # Generate signal
        generatesignal!(replica, replica.isreplica)
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
    noise_val = sum(corr_result[max_idx[1],:])
    PS = pk_val
    N, M = size(corr_result)
    PN = (noise_val - pk_val)/(M - 1)
    # PN = (noise_val - PS)/size(corr_result)[2]
    SNR_est = 10*log10(sqrt(PS/PN))
    return (n0_est, fd_est, SNR_est)
end


#------------------------------------------------------------------------------
#                             OLD IMPLEMENTATON
#------------------------------------------------------------------------------


"""
    courseacquisition(data::GNSSSignal, replica::ReplicaSignal,
                      prn; fd_center=0., fd_range=5000.,
                      fd_rate=0., Δfd=1/replica.t_length,
                      threads=nthreads(), message="Correlating...",
                      operation="replace", start_idx=1,
                      showprogressbar=true)

Initializes `corr_result` array and runs `courseacquisition!` method.
Returns `corr_result` and the course acquired code phase, Doppler
frequency and SNR.
"""
function courseacquisition(data::GNSSSignal, replica::ReplicaSignal,
                           prn; fd_center=0., fd_range=5000.,
                           fd_rate=0., Δfd=1/replica.t_length,
                           threads=nthreads(), message="Correlating...",
                           operation="replace", start_idx=1,
                           showprogressbar=true)
    # Allocate space for correlation result
    corr_result = gencorrresult(fd_range, Δfd, replica.sample_num)
    # Perform course acquisition
    courseacquisition!(corr_result, data, replica, prn;
                       fd_center=fd_center, fd_range=fd_range,
                       fd_rate=fd_rate, Δfd=Δfd, threads=threads)
    n0_est, fd_est, SNR_est = course_acq_est(corr_result)
    return (corr_result, fd_est, n0_est, SNR_est)
end


"""
    courseacquisition!(corr_result::Array{Float64,2},
                       data, replica::ReplicaSignal,
                       prn; fd_center=0., fd_range=5000.,
                       fd_rate=0., Δfd=1/data.t_length,
                       threads=8, message="Correlating...",
                       operation="replace", start_idx=1)

Performs course acquisition on either `GNSSData` or `ReplicaSignal`
type struct using defined `ReplicaSignal` type struct. No need
to use `generatesignal!` before calling this function.
Operates in place on `corr_result`.

*NOTE:* If `data` is a replica signal structure type,
        it's structure type must be the same as the type
        of `replica`. They must also be separately defined
        or one must be a deep copy of the other. Do not
        pass the same structure to the `data` and `replica`
        arguments.

`corr_result` contains |conj(fft(replica)*fft(data)|² per
Doppler bin.
"""
function courseacquisition!(corr_result::Array{Float64,2},
                            data::GNSSSignal, replica::ReplicaSignal,
                            prn; fd_center=0., fd_range=5000.,
                            fd_rate=0., Δfd=1/replica.t_length,
                            message="Correlating...",
                            operation="replace", start_idx=1,
                            showprogressbar=true)
    # Number of data samples
    dsize = replica.sample_num
    # Pre-plan FFTs and IFFTs
    pfft = plan_fft!(replica.data)  # In-place FFT plan
    pifft = plan_ifft!(replica.data) # In-place IFFT plan
    # Carrier wipe data signal, make copy, and take FFT
    datafft = fft(data.data[start_idx:start_idx+dsize-1] .*
                  exp.(-2π.*data.f_if.*data.t[1:dsize].*1im))
    # datafft = fft(data.data)
    # Number of bits representing `data`
    nADC = data.nADC
    # Number of Doppler bins
    doppler_bin_num = Int(fd_range/Δfd*2+1)
    # Loading bar
    if showprogressbar
       p = Progress(doppler_bin_num, 1, message)
    end
    @inbounds for i in 1:doppler_bin_num
        # Calculate Doppler frequency for `i` Doppler bin
        f_d = (fd_center-fd_range) + (i-1)*Δfd
        # Set signal parameters
        definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, ϕ=0., f_if=0.,
                      # include_carrier=true,
                      # include_noise=false,
                      # include_adc=false,
                      code_start_idx=1,
                      nADC=nADC, isreplica=true,
                      noexp=false)
        # Generate signal
        generatesignal!(replica, replica.isreplica)
        # Perform in place FFT on replica
        pfft*replica.data
        # Take conjugate of FFT(replica) and multiply with FFT of
        # data. The result is stored in `replica.signal`
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
        # Update progress bar
        if showprogressbar
          next!(p)
        end
    end
    # Set `isreplica` flag to false in `replica`
    definesignal!(replica, isreplica=false)
    return corr_result
end


"""
    courseacquisition!(corr_result::Array{Float64,2},
                            data::GNSSSignal, replica::ReplicaSignal,
                            prn, N; fd_center=0., fd_range=5000.,
                            fd_rate=0., Δfd=1/replica.t_length,
                            threads=nthreads(), message="Correlating...",
                            operation="add", start_idx=1,
                            showprogressbar=true)

Performs non-coherent course acquisition on either `GNSSData` or
`ReplicaSignal` type struct using defined `ReplicaSignal` type struct.
No need to use `generatesignal!` before calling this function.
Operates in place on `corr_result`.

*NOTE:* If `data` is a replica signal structure type,
        it's structure type must be the same as the type
        of `replica`. They must also be separately defined
        or one must be a deep copy of the other. Do not
        pass the same structure to the `data` and `replica`
        arguments.

`corr_result` contains |conj(fft(replica)*fft(data)|² per
Doppler bin.
"""
function courseacquisition!(corr_result::Array{Float64,2},
                            data::GNSSSignal, replica::ReplicaSignal,
                            prn, N; fd_center=0., fd_range=5000.,
                            fd_rate=0., Δfd=1/replica.t_length,
                            threads=nthreads(), message="Correlating...",
                            operation="add", start_idx=1,
                            showprogressbar=true)
    # Set number of threads to use for FFTW functions
    FFTW.set_num_threads(threads)
    # Number of data samples
    dsize = replica.sample_num
    # Pre-plan FFTs and IFFTs
    pfft = plan_fft!(replica.data)  # In-place FFT plan
    pifft = plan_ifft!(replica.data) # In-place IFFT plan
    # Number of bits representing `data`
    nADC = data.nADC
    # Number of Doppler bins
    doppler_bin_num = Int(fd_range/Δfd*2+1)
    # Carrier wipe data signal, make copy, and take FFT
    datafft = Array{Complex{Float64}}(undef, data.sample_num)
    # Wipeoff f_if from data
    @threads for i in 1:data.sample_num
        @inbounds datafft[i] = data.data[i]*exp(-2π*data.f_if*(data.t[i])*1im)
    end
    # Take in place FFT of each data segment of `replica.t_length` in `datafft`
    @inbounds for n in 1:N
        start_idx = (n-1)*dsize + 1
        pfft*view(datafft, start_idx:start_idx+dsize-1)
    end
    # Create array to hold signal replicas for each Doppler bin
    replicas = Array{Complex{Float64}}(undef, doppler_bin_num, replica.sample_num)
    # Loading bar
    if showprogressbar
       p = Progress(doppler_bin_num*(N+1), 1, message)
    end
    @inbounds for i in 1:doppler_bin_num
        # Calculate Doppler frequency for `i` Doppler bin
        f_d = (fd_center-fd_range) + (i-1)*Δfd
        # Set signal parameters
        definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, ϕ=0., f_if=0.,
                      include_carrier=true,
                      include_noise=false,
                      include_adc=false,
                      code_start_idx=1,
                      nADC=nADC, isreplica=true,
                      noexp=false)
        # Generate signal
        generatesignal!(replica)
        # Perform in place FFT on replica
        pfft*replica.data
        # Take conjugate of replica
        conjA!(replica.data, dsize)
        # Save replica to current Doppler bin
        replicas[i,:] = deepcopy(replica.data)
        # Update progress bar
        if showprogressbar
            next!(p)
        end
    end
    repl = Array{Complex{Float64}}(undef, replica.sample_num)
    @inbounds for n in 1:N
        start_idx = (n-1)*dsize + 1
        datasegment = view(datafft, start_idx:start_idx+dsize-1)
        @inbounds for i in 1:doppler_bin_num
            repl[:] = deepcopy(replicas[i,:])
            # Take product of replica and data
            # Store result in `repl`
            AmultB1D!(repl, datasegment, dsize)
            # Take IFFT in place
            pifft*repl
            # Take `abs2` in place and save result
            # Either replace, add to, or multiply by pre-existing value
            # in `corr_result`
            @threads for j in 1:dsize
                if operation == "replace"
                    @inbounds corr_result[i,j] = abs2(repl[j])
                elseif operation == "add"
                    @inbounds corr_result[i,j] += abs2(repl[j])
                elseif operation == "multiply"
                    @inbounds corr_result[i,j] *= abs2(repl[j])
                end
            end
            # Update progress bar
            if showprogressbar
              next!(p)
            end
        end
    end
    # Set `isreplica` flag to false in `replica`
    definesignal!(replica, isreplica=false)
    return corr_result
end
