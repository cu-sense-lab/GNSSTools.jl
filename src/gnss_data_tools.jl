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
	data, end_idx, dtype = readdatafile!(data, data_type, file_name,
                                         sample_num, start_data_idx)
	# Calculate total data length in seconds
	total_data_length = (end_idx)/f_s
	# Generate time vector
	t = calctvector(sample_num, f_s)
	return GNSSData(file_name, f_s, f_if, t_length, start_data_idx,
                    t, float.(data), String(dtype), data_start_time,
                    site_loc_lla, sample_num, total_data_length, nADC)
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
    return (data, end_idx, :sc8)
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
    return (data, end_idx, :sc4)
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
	data_info_from_name(file_name)

Determine the IF and sampling frequency of the data
from the file name. Looks for the following strings
inside the file name:

- `g1b1`: 25 Msps, centered at 1.57 GHz which contains GPS L1/Galileo E1bc/Beidou B1i
- `g2r2`: 25 Msps, centered at 1.2375 GHz which contains GPS L2 and GLONASS L2
	* `NOTE SUPPORTED`
- `g5`: 25 Msps, centered at 1.17645 GHz which contains GPS L5/Galileo E5a

Returns the IF and sampling frequency in Hz in the order of (f_if, f_s, f_center).
"""
function data_info_from_name(file_name)
	if occursin("g1b1", file_name)
		return g1b1()
	elseif occursin("g2r2", file_name)
		error("L2 band nav signals not supported. Use either L1 (g1b1) or L5 (g5) instead.")
	elseif occursin("g5", file_name)
		return g5()
	else
		error("Cannot determine f_s, f_if, & f_center. Manual specify f_s and f_if.")
	end
end
