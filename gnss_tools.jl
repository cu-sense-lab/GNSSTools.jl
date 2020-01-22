using ProgressMeter
using FFTW
include("signal_generator.jl")


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
struct GNSSData{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,
                A11,A12,A13}
	file_name::A1
	f_s::A2
	f_if::A3
	t_length::A4
	start_data_idx::A5
	t::A6
	data::A7
	data_type::A8
	data_start_time::A9
	site_loc_lla::A10
	sample_num::A11
	total_data_length::A12
	nADC::A13
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
	@inbounds for i in Asize[1]
		@inbounds for j in Asize[2]
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
function AmultB21D!(A, B, Asize=size(A)[1])
	@threads for i in Asize
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
    gencorrresult(fd_range, Δfd, sample_num)

Generates a 2D `Array{Float64,2}` to store the
cross correlation results output from courseacquisition
"""
function gencorrresult(fd_range, Δfd, sample_num)
	return Array{Float64}(undef, Int(2*fd_range/Δfd+1), sample_num)
end


"""
    courseacquisition!(corr_result::Array{Float64,2},
                       data::L5QSignal, replica::L5QSignal,
                       prn; fd_center=0., fd_range=5000.,
                       fd_rate=0., Δfd=1/data.t_length,
                       threads=1)

Performs course acquisition on `L5QSignal` type struct
using defined `L5QSignal` type struct. No need to
use `generatesignal!` before calling this function.
Operates in place on `corr_result`.

Use this version for operating on simulated signals
instead of real data.

`NOTE:` `data` and `signal` must be two different
        `L5QSignal` type structs. Do not use the
        same struct for both argmuents.

`corr_result` contains |conj(fft(replica)*fft(data)|² per
Doppler bin.
"""
function courseacquisition!(corr_result::Array{Float64,2},
                            data::L5QSignal, replica::L5QSignal,
                            prn; fd_center=0., fd_range=5000.,
                            fd_rate=0., Δfd=1/data.t_length,
                            threads=4, message="Correlating...")
	# Set number of threads to use for FFTW functions
	FFTW.set_num_threads(threads)
	# Pre-plan FFTs and IFFTs
	pfft = plan_fft!(replica.signal)  # In-place FFT plan
	pifft = plan_ifft!(replica.signal) # In-place IFFT plan
	# Carrier wipe data signal, make copy, and take FFT
	datafft = fft(data.signal.*exp.(-2π.*data.f_if.*data.t.*1im))
	# Number of bits representing `data`
	nADC = data.nADC
	# Number of data samples
	dsize = data.sample_num
	# Number of Doppler bins
	doppler_bin_num = Int(fd_range/Δfd*2+1)
	# Loading bar
	p = Progress(doppler_bin_num, 1, message)
	@inbounds for i in 1:doppler_bin_num
		# Calculate Doppler frequency for `i` Doppler bin
		f_d = (fd_center-fd_range) + (i-1)*Δfd
		# Set signal parameters
		definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, ϕ=0., f_if=0.,
                      include_carrier=true,
                      include_carrier_amplitude=false,
                      include_noise=false,
                      include_adc=false,
                      nADC=nADC, isreplica=true)
		# Generate signal
		generatesignal!(replica)
		# Perform in place FFT on replica
		pfft*replica.signal
		# Take conjugate of FFT(replica) and multiply with FFT of
		# data. The result is stored in `replica.signal`
		conjAmultB1D!(replica.signal, datafft, dsize)
		# Take IFFT in place and save into `corr_result`
		corr_result[i,:] = abs2.(pifft*replica.signal)
		# Update progress bar
		next!(p)
	end
	# Set `isreplica` flag to false in `replica`
	definesignal!(replica, isreplica=false)
	return corr_result
end


"""
    courseacquisition!(corr_result::Array{Float64,2},
                       data::GNSSData, replica::L5QSignal,
                       prn; fd_center=0., fd_range=5000.,
                       fd_rate=0., Δfd=1/data.t_length,
                       threads=1)

Performs course acquisition on `Data` type struct
using defined `L5QSignal` type struct. No need to
use `generatesignal!` before calling this function.
Operates in place on `corr_result`. 

`corr_result` contains |conj(fft(replica)*fft(data)|² per
Doppler bin.
"""
function courseacquisition!(corr_result::Array{Float64,2},
                            data::GNSSData, replica::L5QSignal,
                            prn; fd_center=0., fd_range=5000.,
                            fd_rate=0., Δfd=1/data.t_length,
                            threads=1, message="Correlating...")
	# Set number of threads to use for FFTW functions
	FFTW.set_num_threads(threads)
	# Pre-plan FFTs and IFFTs
	prfft = fft_plan!(replica.signal)  # In-place FFT plan
	pifft = ifft_plan!(replica.signal) # In-place IFFT plan
	# Carrier wipe data signal, make copy, and take FFT
	datafft = fft(data.data.*exp.(-2π.*data.f_if.*data.t.*1im))
	# Number of bits representing `data`
	nADC = data.nADC
	# Number of data samples
	dsize = data.sample_num
	# Number of Doppler bins
	doppler_bin_num = Int(fd_range/Δfd*2+1)
	# Loading bar
	p = Progress(doppler_bin_num, 1, message)
	@inbounds for i in 1:doppler_bin_num
		println(i)
		# Calculate Doppler frequency for `i` Doppler bin
		f_d = (fd_center-fd_range) + (i-1)*Δfd
		# Set signal parameters
		definesignal!(replica;
                      prn=prn, f_d=f_d,
                      fd_rate=fd_rate, ϕ=0., f_if=0.,
                      include_carrier=true,
                      include_carrier_amplitude=false,
                      include_noise=false,
                      include_adc=false,
                      nADC=nADC)
		# Generate signal
		generatesignal!(replica)
		# Perform in place FFT on replica
		prfft*replica.signal
		# Take conjugate of FFT(replica) and multiply with FFT of
		# data. The result is stored in `replica.signal`
		conjAmultB1D!(replica.signal, datafft, dsize)
		# Take IFFT in place and save into `corr_result`
		corr_result[i,:] = abs2.(pifft*replica.signal)
		# Update progress bar
		next!(p)
	end
	return corr_result
end


