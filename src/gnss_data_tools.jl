"""
    loaddata(data_type, file_name, f_s, f_center, f_gnss, t_length;
             start_data_idx=1, site_lla=missing, data_start_time=missing)


Loads data from `scN` file and loads into `GNSSData` type struct.

`site_lla` is `[latitude, longitude, height]` in (degrees, degrees, meters).

`data_start_time` is `[year, month, day, hour, minute, seconds]`, where "ss" is
decimal seconds.


Required Arguments:

- `data_type`: one of the following for `Complex{IntN}` datasets
 	* `:sc4`
	* `:sc8`
	* `:sc16`
	* `:sc32`
	* `:sc64`
- `file_name::String`: data file name
- `f_s::Float64`: sampling rate of data in Hz
- `f_center::Float64`: center frequency of the receiver in Hz
- `f_gnss::Float64`: center frequency of the GNSS signal in Hz
- `t_length::Float64`: the length of data in seconds to be loaded


Optional Arguments:

- `skip_to::Float64`: where the first data sample will be loaded in the data
                      `(default = 0.1sec)`
	* units are in seconds
- `site_lla`: receiver LLA location in the format `(lat, lon, height)`
	* `lat`, latitude, is in degrees
	* `lon`, longitude, is in degrees
	* `height` is in meters
	* `(default = missing)`
- `data_start_time`: data UTC start time `(year, month, day, hour, minute, seconds)`
	* `seconds` can be `Int` or `Float64`, all others are `Int`
	* `(default = missing)`


Returns:

- `GNSSData` struct
"""
function loaddata(data_type, file_name, f_s, f_center, f_gnss, t_length;
                  skip_to=0.1, site_loc_lla=missing,
                  data_start_time=missing)
	# Compute number of samples to extract
	sample_num = Int(f_s * t_length)
	# Calculate f_if
	f_if = f_gnss - f_center
	# Determine nADC
	if data_type == :sc4
		data_type = Val(:sc4)
		nADC = 4
		data = Array{Complex{Int8}}(undef, sample_num)
	elseif data_type == :sc8
		data_type = Val(:sc8)
		nADC = 8
		data = Array{Complex{Int8}}(undef, sample_num)
	elseif data_type == :sc16
		data_type = Val(:sc16)
		nADC = 16
		data = Array{Complex{Int16}}(undef, sample_num)
	elseif data_type == :sc32
		data_type = Val(:sc32)
		nADC = 32
		data = Array{Complex{Int32}}(undef, sample_num)
	elseif data_type == :sc64
		data_type = Val(:sc64)
		nADC = 64
		data = Array{Complex{Int64}}(undef, sample_num)
	else
		error("Invalid data type: $(typeof(data_type))")
	end
	# Read data
	start_data_idx = floor(Int, f_s*skip_to)
	data, end_idx, dtype = readdatafile!(data, data_type, file_name,
	                                     sample_num, nADC, start_data_idx)
	# Calculate total data length in seconds
	total_data_length = (end_idx)/f_s
	# Generate time vector
	t = calctvector(sample_num, f_s)
	return GNSSData(file_name, f_s, f_if, t_length, start_data_idx,
	                t, float.(data), String(dtype), data_start_time,
	                site_loc_lla, sample_num, total_data_length, nADC)
	# return GNSSData(file_name, f_s, f_if, t_length, start_data_idx,
	#                 t, data, String(dtype), data_start_time,
	#                 site_loc_lla, sample_num, total_data_length, nADC)
end


"""
	reloaddata!(gnss_data::GNSSData, start_data_idx, sample_num)


Reloads portion of data up to the length of the data array inside `gnss_data`.


Required Arguments:

- `gnss_data::GNSSData`: structure that holds loaded data from data files
- `start_data_idx::Int`: where the first data sample will be loaded in the data
- `sample_num::Int`: number of samples after `start_data_idx` to load into
                     `gnss_data`
    * must not be greater than `gnss_data.sample_num`


Modifies in Place and Returns:

- `gnss_data::GNSSData`
"""
function reloaddata!(gnss_data::GNSSData, start_data_idx,
	                 nADC, sample_num=gnss_data.sample_num)
	file_name = gnss_data.file_name
	data_type = gnss_data.data_type
	data = gnss_data.data
	if sample_num <= gnss_data.sample_num
		data, end_idx, dtype = readdatafile!(data, Val(Symbol(data_type)),
		                                     file_name, sample_num,
											 start_data_idx)
        f = open(file_name, "r")
		byte_num = Int(2*nADC/8)
		seek(f, byte_num*start_data_idx)
		read!(f, data)
		close(f)
	else
		error("`sample_num` is greater than data sample size!")
	end
	return gnss_data
end


"""
	readdatafile!(data, data_type, file_name, sample_num, start_idx=1)


Loads `scN` data files. First N-bit number is real, second is imaginary.


Required Arguments:

- `data::Vector{Complex{IntN}}`: the complex data vector to store the new data
- `data_type::Val{:scN}`: only accepts `Val(:scN)` to indicate data is N-bit
                          complex
- `file_name::String`: data file name
- `sample_num::Int`: `[NOT USED]` number of samples after `start_data_idx` to
                     load


Optional Arguments:

- `start_idx::Int`: where the first data sample will be loaded in the data
                    `(default = 0)`


Returns:

- `data::Vector{Complex{IntN}}`: elements are replaced with new data
- `end_idx::Int`: the last sample index loaded plus one
- `:scN::Symbol`: symbol which describes the type of data
"""
function readdatafile!(data, data_type, file_name, sample_num,
                       nADC, start_idx=0)
	# Open file
	f = open(file_name, "r")
	# Go to start location
	byte_num = 2*Int(nADC/8)
	seek(f, byte_num*start_idx)
	# Load Complex{IntN} values directly from file
	read!(f, data)
	# Get the index value for the end of the file
	end_idx = Int(position(seekend(f))/byte_num)
    close(f)
    return (data, end_idx, gnsstypes[data_type])
end


"""
	readdatafile!(data, data_type::Val{:sc4}, file_name, sample_num,
				  start_idx=1)


Loads `sc4` data files. For each UInt8 number, the LSB is real and MSB is
imaginary.

**NOTE:** this method is slower than the version used for 8-bit complex data.


Required Arguments:

- `data::Vector{Complex{Int}}`: the complex data vector to store the new data
- `data_type::Val{:sc4}`: only accepts `Val(:sc4)` to indicate data is 4-bit
                          complex
- `file_name::String`: data file name
- `sample_num::Int`: `[NOT USED]` number of samples after `start_data_idx` to
                     load


Optional Arguments:

- `start_idx::Int`: where the first data sample will be loaded in the data
	                `(default = 0)`


Returns:

- `data::Vector{Complex{Int}}`: elements are replaced with new data
- `end_idx::Int`: the last sample index loaded plus one
- `:sc4::Symbol`: symbol which describes the type of data
"""
function readdatafile!(data, data_type::Val{:sc4}, file_name, sample_num,
                       nADC, start_idx=0)
	# Open file
	f = open(file_name, "r")
	# Go to start location
	seek(f, start_idx)
	# Load UInt8 values and convert to Complex{Int8}
	data[:] = bytetocomplex.(read(f, sample_num))
	# Get the index value for the end of the file
	end_idx = position(seekend(f))
    close(f)
    return (data, end_idx, :sc4)
end


"""
    bytetocomplex(byte::UInt8)


Converts a `UInt8` value to `Complex{Int8}`. Assumes that the LSB and MSB are the
real and imaginary parts, respectively.


Required Arguments:

- `byte::UInt8`: a single byte to converted to `Complex{Int8}`


Returns:

- `Complex{Int8}`: complex number converted from `byte`
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
	FileInfo


`Struct` that holds data file information based off its file name.


Fields:

- `f_s::Float64`: sampling rate of signal in Hz
- `f_if::Float64`: signal IF frequency in Hz
- `f_center::Float64`: receiver center frequency in Hz
- `sig_freq::Float64`: navigation band frequency in Hz
                       (e.g. for L1, `sig_freq = L1_freq`)
- `sigtype::T1`: `[DEPRICATED]`
- `data_type::T2`: the data type of the data file which is either set to
                   `Val(:sc4)` or `Val(:sc8)`
- `timestamp::T3`: six element `Tuple` with format
                   `(year, month, day, hour, minute, second)` where everything
				   other than `second`, which can be either `Int` or `Float64`,
				   is an Int
- `timestamp_JD::Float64`: the Julia date of the timestamp
- `file_name::String`: file name of the data file
"""
struct FileInfo{T1,T2,T3}
	f_s::Float64
	f_if::Float64
	f_center::Float64
	sig_freq::Float64
	sigtype::T1
	data_type::T2
	timestamp::T3
	timestamp_JD::Float64
	file_name::String
end


"""
	data_info_from_name(file_name)


Determine the IF and sampling frequency of the data from the file name. Looks
for the following strings inside the file name:

- `g1b1`: 25 Msps, centered at 1.57 GHz which may contain the followin signals
	* GPS L1
	* Galileo E1bc
	* Beidou B1i
- `g2r2`: `[NOTE SUPPORTED]` 25 Msps, centered at 1.2375 GHz which may contain
          the followin signals
	* GPS L2
	* GLONASS L2
- `g5`: 25 Msps, centered at 1.17645 GHz which may contain the followin signals
	* GPS L5
	* Galileo E5a

Also looks at the file name extension to determine the type of data stored in
the file where

- `.sc4`: 4-bit complex data where for each UInt8 number in the data file, the
  LSB is real and MSB is imaginary
- `.sc8`: 8-bit complex data where for each pair of Int8 numbers, the first
  value is real and the second is imaginary


Required Arguments:

- `file_name::String`: the data file name to parse


Returns:

- `FileInfo` struct
"""
function data_info_from_name(file_name)
	# Get file type (8- or 4-bit complex)
	if occursin("sc8", file_name)
		data_type = Val(:sc8)
	elseif occursin("sc4", file_name)
		data_type = Val(:sc4)
	else
		error("File type not supported. Only sc8 and sc4 files are supported.")
	end
	# Determine sampling and IF frequency and frequency center
	f_s, f_if, f_center, sig_freq, sigtype = get_signal_type(file_name)
	timestamp, timestamp_JD = find_and_get_timestamp(file_name)
	return FileInfo(f_s, f_if, f_center, sig_freq, sigtype, data_type,
	                timestamp, timestamp_JD, file_name)
end
