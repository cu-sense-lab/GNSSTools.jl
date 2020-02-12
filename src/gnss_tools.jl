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

