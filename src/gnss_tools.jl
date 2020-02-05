"""
    fft_correlate(data, reference)

Calculate the cyclical FFT based correlation
between the data and the reference signal.

Returns:

- Array containing the correlation result
"""
function fft_correlate(data, reference)
    return ifft(conj!(fft(reference)).*fft(data))
end


"""
    Data()


Structure for holding signal data.
"""
struct GNSSData{A1,A2,A3}
	file_name::String
	f_s::Float64
	f_if::Float64
	t_length::Float64
	start_data_idx::Int64
	t::Array{Float64,1}
	data::Array{Complex{Float64},1}
	data_type::A1
	data_start_time::A2
	site_loc_lla::A3
	sample_num::Int64
	total_data_length::Float64
	nADC::Int64
end


"""
    loaddata(data_type, file_name, f_s, f_if, t_length;
             start_data_idx=1, site_lla=missing, data_start_time=missing)

Loads data from sc8 file and loads into GNSSData type struct.

site_lla is [latitude, longitude, height] in degrees and meters

data_start_time is [YYYY, MM, DD, HH, mm, ss], where "ss" is decimal
seconds.
"""
function loaddata(data_type, file_name, f_s, f_if, t_length;
                  start_data_idx=1, site_loc_lla=missing,
                  data_start_time=missing)
	# Determine nADC
	if typeof(data_type) == Val{:sc8}
		nADC = 8
	elseif typeof(data_type) == Val{:sc4}
		nADC = 4
	else
		error("Invalid data type: $(typeof(data_type))")
	end
	# Compute number of samples to extract
	sample_num = Int(f_s * t_length)
	# Read data
	data = Array{Complex{Int8}}(undef, sample_num)
	data, end_idx = readdatafile!(data, data_type, file_name,
                                  sample_num, start_data_idx)
	# Calculate total data length in seconds
	total_data_length = (end_idx)/f_s
	# Generate time vector
	t = Array(0:1/f_s:t_length-1/f_s)  # s
	return GNSSData(file_name, f_s, f_if, t_length, start_data_idx,
                    t, float.(data), data_type, data_start_time, site_loc_lla,
                    sample_num, total_data_length, nADC)
end


"""
    readdatafile(data_type::Val{:sc8}, file_name, sample_num,
                 start_idx=0, message="Loading data...")

Loads sc8 data files. First 8-bit number is real, second is
imaginary.
"""
function readdatafile!(data, data_type::Val{:sc8}, file_name, sample_num,
                       start_idx=0, message="Loading data...")
	# Open file
	f = open(file_name, "r")
	# Go to start location
	seek(f, 2*(start_idx-1))
	# Load Complex{Int8} values directly from file
	read!(f, data)
	# Get the index value for the end of the file
	end_idx = position(seekend(f)) + 1
    close(f)
    return (data, end_idx)
end


"""
    readdatafile(data_type::Val{:sc4}, file_name, sample_num,
                 start_idx=0, message="Loading data...")

Loads sc8 data files. For each UInt8 number,
the LSB is real and MSB is imaginary.
"""
function readdatafile!(data, data_type::Val{:sc4}, file_name, sample_num,
                       start_idx=0, message="Loading data...")
	# Open file
	f = open(file_name, "r")
	# Go to start location
	seek(f, start_idx-1)
	# Load UInt8 values and convert to Complex{Int8}
	data[:] = bytetocomplex.(read(f, sample_num))
	# Get the index value for the end of the file
	end_idx = position(seekend(f)) + 1
    close(f)
    return (data, end_idx)
end


"""
    bytetocomplex(byte::UInt8)

Converts a UInt8 value to Complex{Int8}
"""
function bytetocomplex(byte::UInt8)
	# Get unsigned I and Q values
	Iuns = byte & 0x0f
	Quns = byte >> 0x04
	# Get I and Q signs
	Isign = (Iuns&0x08) >> 0x03
	Qsign = (Quns&0x08) >> 0x03
	# Get I and Q magnitudes
	Imag = xor(Isign*0x07, Iuns&0x07)
	Qmag = xor(Qsign*0x07, Quns&0x07)
	# Compute the complex number
	return (-1)^(Isign)*Imag + (-1)^(Qsign)*Qmag*1im
end


"""
    AmultB2D!(A, B, Asize=size(A))

Multiply contents of A in place with contents of B.
Both A and B should be 2D arrays and be the same size.
"""
function AmultB2D!(A, B, Asize=size(A))
	@inbounds for i in 1:Asize[1]
		@inbounds for j in 1:Asize[2]
			A[i,j] = A[i,j] * B[i,j]
		end
	end
	return A
end


"""
    AmultB1D!(A, B, Asize=size(A))

Multiply contents of A in place with contents of B.
Both A and B should be 1D arrays and be the same size.
"""
function AmultB1D!(A, B, Asize=size(A)[1])
	@threads for i in 1:Asize
		@inbounds A[i] = A[i] * B[i]
	end
	return A
end


"""
    conjAmultB1D!(A, B, Asize=size(A))

Multiply contents of conj(A) in place with contents of B.
Both A and B should be 1D arrays and be the same size.
"""
function conjAmultB1D!(A, B, Asize=size(A)[1])
	@threads for i in 1:Asize
		@inbounds A[i] = conj(A[i]) * B[i]
	end
	return A
end


"""
    conjA!(A, Asize=size(A))

Takes the conjugate of A in place.
A should be a 1D array.
"""
function conjA!(A, Asize=size(A)[1])
    @threads for i in 1:Asize
        @inbounds A[i] = conj(A[i])
    end
    return A
end


"""
    gencorrresult(fd_range, Δfd, sample_num)

Generates a 2D `Array{Float64,2}` to store the
cross correlation results output from courseacquisition
"""
function gencorrresult(fd_range, Δfd, sample_num; iszeros=false)
    if iszeros
        return zeros(Float64, Int(2*fd_range/Δfd+1), sample_num)
    else
	   return Array{Float64}(undef, Int(2*fd_range/Δfd+1), sample_num)
    end
end


"""
    courseacquisition!(corr_result::Array{Float64,2},
                       data, replica::L5QSignal,
                       prn; fd_center=0., fd_range=5000.,
                       fd_rate=0., Δfd=1/data.t_length,
                       threads=8, message="Correlating...",
                       operation="replace", start_idx=1)

Performs course acquisition on either `GNSSData` or `L5QSignal`
type struct using defined `L5QSignal` type struct. No need
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
                            data, replica::L5QSignal,
                            prn; fd_center=0., fd_range=5000.,
                            fd_rate=0., Δfd=1/replica.t_length,
                            threads=8, message="Correlating...",
                            operation="replace", start_idx=1,
                            showprogressbar=true)
	# Set number of threads to use for FFTW functions
	FFTW.set_num_threads(threads)
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
                      include_carrier=true,
                      include_noise=false,
                      include_adc=false,
                      nADC=nADC, isreplica=true)
		# Generate signal
		generatesignal!(replica)
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
                       data, replica::L5QSignal,
                       prn, N=1; fd_center=0., fd_range=5000.,
                       fd_rate=0., Δfd=1/data.t_length,
                       threads=8, message="Correlating...",
                       operation="replace", start_idx=1)

Performs non-coherent course acquisition on either `GNSSData` or
`L5QSignal` type struct using defined `L5QSignal` type struct.
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
                            data, replica::L5QSignal,
                            prn, N; fd_center=0., fd_range=5000.,
                            fd_rate=0., Δfd=1/replica.t_length,
                            threads=8, message="Correlating...",
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
    # Generate all replicas and store in `replicas`
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
                      nADC=nADC, isreplica=true)
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


"""
    calcsnr(x)

Calculates the SNR of the correlation peak in `x`.
"""
function calcsnr(x)
	N = length(x)
	amplitude = sqrt(maximum(abs2.(x)))
	PS = 2*amplitude^2
	PN = 0.
	@threads for i in 1:N
		@inbounds PN += abs2(x[i])
	end
	PN -= PS/(N-2)
	return 10*log10(PS/PN)
end

