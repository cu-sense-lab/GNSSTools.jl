"""
    AmultB2D!(A, B, Asize=size(A))


Multiply contents of A in place with contents of B. Both A and B should be 2D
arrays and be the same size.


Required Arguments:

- `A::Array{T,2}`: 2D dimmensional array, `A`, with element type `T`
- `B::Array{T,2}`: 2D dimmensional array, `B`, with element type `T`


Optional Arguments:

- `Asize::Tuple(Int,Int)`: size of array `A`, assuming that `B` is the same size
                           size as `A` `(default = size(A))`


Modifies in Place and Returns:

- `A::Array{T,2}`
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
    AmultB1D!(A::Vector, B::Vector, Asize=size(A)[1])


Multiply contents of A in place with contents of B. Both A and B should be 1D
arrays and be the same size.


Required Arguments:

- `A::Array{T,1}`: 1D dimmensional array, `A`, with element type `T`
- `B::Array{T,1}`: 1D dimmensional array, `B`, with element type `T`


Optional Arguments:

- `Asize::Int`: size of array `A`, assuming that `B` is the same size
                size as `A` `(default = size(A)[1])`


Modifies in Place and Returns:

- `A::Array{T,1}`
"""
function AmultB1D!(A::Vector, B::Vector, Asize=size(A)[1])
    @threads for i in 1:Asize
        @inbounds A[i] = A[i] * B[i]
    end
    return A
end


"""
    conjAmultB1D!(A, B, Asize=size(A))


Multiply contents of conj(A) in place with contents of B. Both A and B should
be 1D arrays and be the same size.


Required Arguments:

- `A::Array{Complex{T},1}`: 1D dimmensional array, `A`, with element type `T`
- `B::Array{Complex{T},1}`: 1D dimmensional array, `B`, with element type `T`


Optional Arguments:

- `Asize::Int`: size of array `A`, assuming that `B` is the same size
                size as `A` `(default = size(A)[1])`


Modifies in Place and Returns:

- `A::Array{Complex{T},1}`
"""
function conjAmultB1D!(A::Vector, B::Vector, Asize=size(A)[1])
    @threads for i in 1:Asize
        @inbounds A[i] = conj(A[i]) * B[i]
    end
    return A
end


"""
    conjA!(A, Asize=size(A))


Takes the conjugate of A in place. A should be a 1D array.


Required Arguments:

- `A::Array{Complex{T},1}`: 1D dimmensional array, `A`, with element type `T`


Optional Arguments:

- `Asize::Int`: size of array `A` `(default = size(A)[1])`


Modifies in Place and Returns:

- `A::Array{Complex{T},1}`
"""
function conjA!(A::Vector, Asize=size(A)[1])
    @threads for i in 1:Asize
        @inbounds A[i] = conj(A[i])
    end
    return A
end


"""
    calcsnr(x)


Calculates the SNR of the correlation peak in `x`.


Required Arguments:

- `x::Array{T,1}`: 1 dimmenstional array of type `T`


Returns:

- `Float64`: the SNR of the peak in `x` in dB
"""
function calcsnr(x::Vector)
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


"""
    fft_correlate(data, reference)


Calculate the cyclical FFT based correlation between the data and the reference
signal.


Required Arguments:

- `data::Array{T,N}`: original `N` dimmensional data signal of element
                      type `T`
- `reference::Array{T,N}`: `N` dimmensional reference signal of element
                           type `T`


Returns:

- `Array{Complex{Float64},N}`: `N` dimmensional array containing the correlation
                               result
"""
function fft_correlate(data, reference)
    return ifft(conj!(fft(reference)).*fft(data))
end


"""
    gnsstypes


Dictionary containing the qyuivalent strings for each type used in `GNSSTools`.
"""
const gnsstypes = Dict(Val{:l5q}() => "l5q",
                       Val{:l5i}() => "l5i",
                       Val{:l1ca}() => "l1ca",
                       Val{:fft}() => "fft",
                       Val{:carrier}() => "carrier",
					   Val{:sc4}() => "sc4",
					   Val{:sc8}() => "sc8",
					   Val{:sc16}() => "sc16",
                       Val{:sc32}() => "sc32",
					   Val{:sc64}() => "sc64")


"""
	calcinitcodephase(code_length::Int, f_code_d::Number, f_code_dd::Number,
					  f_s::Number, code_start_idx::Number)


Calculates the initial code phase of a given code where f_d and fd_rate are the
Doppler affected code frequency and code frequency rate, respectively. The
returned value is the fractional index location in the original code vector in
`CodeType.codes[i]`, where `i` is the current code being processed.


Required Arguments:

**NOTE:** Arguments can be either `Float64` or `Int`.

- `code_length::Int`: number of bits in code
- `f_code_d::Number`: Doppler adjusted code chipping rate in Hz
- `f_code_dd::Number`: Doppler rate adjusted code chipping rate rate in Hz
- `f_s::Number`: sampling rate of data in Hz
- `code_start_idx::Number`: index in `ReplicaSignal.data` where all layers of
                            code start


Returns:

- `Float64`: Fractional code index location
"""
function calcinitcodephase(code_length::Int, f_code_d::Number,
	                       f_code_dd::Number, f_s::Number,
						   code_start_idx::Number)
    t₀ = (code_start_idx-1)/f_s
    init_phase = -f_code_d*t₀ - 0.5*f_code_dd*t₀^2
    return (init_phase%code_length + code_length)%code_length
end


"""
	calccodeidx(init_chip::Number, f_code_d::Number, f_code_dd::Number,
				t::Number, code_length::Int)


Calculates the index in the codes for a given t.


Required Arguments:

- `init_chip::Number`: initial code phase at `t=0s` returned from
                       `calcinitcodephase`
- `f_code_d::Number`: Doppler adjusted code chipping rate in Hz
- `f_code_dd::Number`: Doppler rate adjusted code chipping rate rate in Hz
- `t::Number`: current time elapsed in seconds
- `code_length::Int`: number of bits in code


Returns:

- `Int`: the current code index at a given `t`
"""
function calccodeidx(init_chip::Number, f_code_d::Number, f_code_dd::Number,
	                 t::Number, code_length::Int)
    return Int(floor(init_chip+f_code_d*t+0.5*f_code_dd*t^2)%code_length)+1
end


"""
	calc_doppler_code_rate(f_code::Number, f_carrier::Number, f_d::Number,
						   fd_rate::Number)


Calculates the adjusted code chipping rate and chipping rate rate based off the
Doppler and Doppler rate.


Required Arguments:

- `f_code::Number`: chipping rate of code with no Doppler and Doppler rate in Hz
- `f_carrier::Number`: carrier frequency of signal in Hz
- `f_d::Number`: Doppler frequency in Hz
- `fd_rate::Number`: Doppler frequency rate in Hz


Returns:

- `(f_code_d::Float64, f_code_dd::Float64)::Tuple` where
	* `f_code_d`: Doppler adjusted code chipping rate in Hz
	* `f_code_dd`: Doppler rate adjusted code chipping rate rate in Hz
"""
function calc_doppler_code_rate(f_code::Number, f_carrier::Number, f_d::Number,
	                            fd_rate::Number)
	f_code_d = f_code*(1. + f_d/f_carrier)
	f_code_dd = f_code*fd_rate/f_carrier
	return (f_code_d, f_code_dd)
end


"""
    calctvector(N::Int, f_s::Number)


Generates an `N` long time vector with time spacing `Δt` or `1/f_s`.


Required Arguments:

- `N::Int`: length of the time vector
- `f_s::Number`: sampling rate of signal in Hz


Returns:

- `t::Array{Float64,1}`: time vector of length `N` and spacing of `1/f_s`
                         seconds
"""
function calctvector(N::Int, f_s::Number)
    # Generate time vector
    t = Array{Float64}(undef, N)
    @threads for i in 1:N
        @inbounds t[i] = (i-1)/f_s
    end
    return t
end


"""
    meshgrid(x::Vector, y::Vector)


Generate a meshgrid the way Numpy would in Python.


Required Arguments:

- `x::Vector{T}`, 1D array of element type `T`
- `y::Vector{T}`, 1D array of element type `T`


Returns:

- `(X, Y)::Tuple` where
	* `X::Array{eltype(x),2}`: 2D array of size `(length(y), length(x))`
	* `Y::Array{eltype(y),2}`: 2D array of size `(length(y), length(x))`
"""
function meshgrid(x::Vector, y::Vector)
    xsize = length(x)
    ysize = length(y)
    X = Array{eltype(x)}(undef, ysize, xsize)
    Y = Array{eltype(y)}(undef, ysize, xsize)
    for i in 1:ysize
        for j in 1:xsize
            X[i,j] = x[j]
            Y[i,j] = y[i]
        end
    end
    return (X, Y)
end


"""
	find_sequence(file_name, digit_nums, separators=missing)


Find a number sequence in a file name. `digit_nums` can be a number or array.
`separators` can also be either a number or array containing the characters such
as `_` and `.` that may be between the numbers. Returns the string containing
the sequence if found. Returns `missing` if not found.


Required Arguments:

- `file_name`: name of file
- `digit_nums`: can be an array or single number indicating number of digits in
                `file_name` to find
    * e.g. to find `1970-01-01`, user would pass `digit_nums` as `[4, 2, 2]`


Optional Arguments:

- `separators`: the characer used to separate the number sequence
	* e.g. to find `1970-01-01`, user would pass `separators` as `"-"`
	* `(default = missing)`


Returns:

- `String`: if sequence is found, a copy of the sequence in `file_name` is
            returned
    * returns `missing` if sequence is not found
"""
function find_sequence(file_name, digit_nums, separators=missing)
	total_digits = sum(digit_nums)
	section_num = length(digit_nums)
	current_section = 1
	sequence_found = false
	sequence_complete = false
	between_sections = false
    sequence_counter = 0
	sequence_start = 1
    sequence_stop = 1
	for i in 1:length(file_name)
        if isdigit(file_name[i])
            if sequence_found == false
                sequence_start = i
            end
            sequence_found = true
            sequence_counter += 1
			if between_sections
				current_section += 1
				between_sections = false
			end
        else
            if sequence_found && (sequence_counter == sum(digit_nums[1:current_section])) && occursin(file_name[i], separators)
				between_sections = true
            else
                sequence_found = false
                sequence_counter = 0
                sequence_idx_start = 1
                sequence_idx_stop = 1
				current_section = 1
				between_sections = false
            end
        end
		if (sequence_counter == total_digits) && (current_section == section_num)
			sequence_stop = i
			sequence_complete = true
			break
		end
    end
	if sequence_complete
		return file_name[sequence_start:sequence_stop]
	else
		return missing
	end
end


"""
    find_and_get_timestamp(file_name)


Find a sequency of 8 digits with `_` separating it from a sequence of 6 digits.
Return the timestamp tuple. Expected sequence to find is `YYYYMMDD_HHMMSS`.


Required Arguments:

- `file_name`: name of data file


Returns:

- `timestamp`: `Tuple` containing date and time
	* format is `(year, month, day, hour, minute, second)`
- `timestamp_JD`: Julian date of `timestamp`
"""
function find_and_get_timestamp(file_name)
	timestamp_string = find_sequence(file_name, [8, 6], "_")
    if ~ismissing(timestamp_string)
        year = parse(Int, timestamp_string[1:4])
    	month = parse(Int, timestamp_string[5:6])
    	day = parse(Int, timestamp_string[7:8])
    	hour = parse(Int, timestamp_string[10:11])
    	minute = parse(Int, timestamp_string[12:13])
    	second = parse(Int, timestamp_string[14:15])
    	timestamp = (year, month, day, hour, minute, second)
    	timestamp_JD = DatetoJD(timestamp...)
        return (timestamp, timestamp_JD)
    else
        error("Data timestamp not found in file name. Please supply it manually.")
    end
end


"""
	random_databits(chipping_rate, t_length; prns=missing)


Generates single or multiple random bit streams. Multiple bit streams are stored
in a dictionary wit the keys being the contents of `prns`.


Required Arguments:

- `chipping_rate`: databit chipping rate in Hz
- `t_length`: length of databit sequence in seconds


Optioonal Arguments:

- `prns`: vector containing PRN numbers `(default = missing)`


Returns:

- Vector of databits or dictionary of `length(prns)` databit sequences
"""
function random_databits(chipping_rate, t_length; prns=missing)
	if ismissing(prns)
		return rand(0:1, ceil(Int, chipping_rate*t_length))
	else
		databits = Dict{eltype(prns),Vector{Int}}()
		for prn in prns
			databits[prn] = rand(0:1, ceil(Int, chipping_rate*t_length))
		end
		return databits
	end
end


"""
	copy_dictionary(dict, key_list)


Copies a dictionary keys and values with keys defined as `b`.


Required Arguments:

- `dict`: dictionary to copy
- `key_list`: list of keys to copy


Return:

- dictionary containing partial to full copy of key_list
"""
function copy_dictionary(dict, key_list)
	dict_copy = Dict{eltype(keys(dict)),Vector{eltype(dict[key_list[1]])}}()
	for key in key_list
		dict_copy[key] = dict[key]
	end
	return dict_copy
end



"""
    get_signal_type(file_name)


Find the signal type of the data based off its file name. Only checks if signal
type is L1 or L5 signal. L2 is not supported.


Required Arguments:

- `file_name`: name of data file


Returns:

- `f_s`: sampling rate in Hz
- `f_if`: IF frequency in Hz
- `f_center`: receiver center frequency in Hz
- `sig_freq`: signal carrier frequency in Hz
- `sigtype`: either `Val(:l1ca)` or `Val(:l5q)`
"""
function get_signal_type(file_name)
    # Determine sampling and IF frequency and frequency center
	if occursin("g1b1", file_name)
		f_s, f_if, f_center, sig_freq, sigtype = g1b1()
	elseif occursin("g2r2", file_name)
		error("L2 band signals not supported. Use either L1 (g1b1) or L5 (g5) files instead.")
	elseif occursin("g5", file_name)
		f_s, f_if, f_center, sig_freq, sigtype = g5()
	else
		# Must check for information on center frequency and sampling rate
		# in file name. If not present, user must specify manually.
		data_freq_string =  find_sequence(file_name, [4,2,2], "._M")
		if ismissing(data_freq_string)
			data_freq_string =  find_sequence(file_name, [4,2,1], "._M")
		end
		if ismissing(data_freq_string)
			error("Cannot determine f_s, f_if, & f_center. Manually specify f_s and f_if.")
		else
			data_freq_string = replace(data_freq_string, "M"=>"")
			data_freq_string = split(data_freq_string, "_")
			f_center = parse(Float64, data_freq_string[1]) * 1e6  # Hz
			f_s = parse(Float64, data_freq_string[2]) * 1e6  # Hz
			# Determine the smallest difference between `f_center` and
			# all GPS carrier frequencies. The closest difference determines
			# what the signal type is.
			Δf = abs.(f_center .- [L1_freq, L2_freq, L5_freq])
			idx = argmin(Δf)
			if idx == 1
				sigtype = Val(:l1ca)
				sig_freq = L1_freq
			elseif idx == 3
				sigtype = Val(:l5q)
				sig_freq = L5_freq
			else
				error("L2 signals are not supported. Supported signals are L1 and L5.")
			end
			f_if = abs(sig_freq-f_center)
		end
	end
	return (f_s, f_if, f_center, sig_freq, sigtype)
end


"""
	make_subplot(fig, row, col, i; projection3d=false, aspect="auto")


Makes a subplot using the PyPlot module. Provides easier setup for 3D subplots.


Required Arguments:

- `fig`: PyPlot figure object
- `row::Int`:specifies the number of rows in the figure
- `col::Int`:specifies the number of columns in the figure
- `i::Int`: specifies the subplot placement in the figure `(range = (1,row*col))`


Optional Arguments:

- `projection3d`: set to `true` if generating a 3D plot `(default = `false`)`
- `aspect`: subplot aspect ratio; refer to PyPlot `(default = "auto")`
	* does not work if `projection3d` is set to `true`


Returns:

- `fig`: PyPlot figure object
- `ax`: PyPlot subplot axis object
"""
function make_subplot(fig, row, col, i; projection3d=false, aspect="auto")
	if projection3d
		# Must be called in local scope in order to make 3D subplots
		mplot3d = PyPlot.PyObject(PyPlot.axes3D)
		ax = fig.add_subplot(row, col, i, projection="3d")
	else
		ax = fig.add_subplot(row, col, i, aspect=aspect)
	end
	return (fig, ax)
end


"""
    plot_spectrum(x::Vector, f_s, log_freq::Bool=false)


Plot a spectrum of the time domain data `x` with frequency range determined
by the sampling rate, `f_s`.


Required Arguments:

- `x::Vector{T}`: 1D array with element types `T`
- `f_s::Number`: sampling rate of `x`


Optional Arguments:

- `log_freq::Bool`: flag for if plot shows linear frequency or `log10(freqs)`
                    where `freqs` is the array of frequencies spanning
					`(0, f_s/2) Hz` `(default = false)`


Returns `nothing`

Displays plot of frequencies against power of `fft(x)`
"""
function plot_spectrum(x::Vector, f_s::Number; log_freq::Bool=true, 
	                   new_fig=true, log_X=false)
    N = length(x)
    t_length = N/f_s
	X = fft(x)[1:floor(Int, N/2)+1]
	X = (1/(f_s*N)) .* abs2.(X)
	X[2:end-1] .= 2 .* X[2:end-1]
	freqs = Array(range(0, f_s/2, length=length(X)))
	if new_fig
    	fig = figure()
		ax = fig.add_subplot(1, 1, 1)
	end
    if log_freq
		if new_fig
			ax.set_xscale("log")
			if log_X
				ax.set_yscale("log")
			end
		else
			xscale("log")
			if log_X
				yscale("log")
			end
		end
	end
	if log_X
		plot(freqs, X)
	else
		plot(freqs, 10 .* log10.(X))
	end
	xlabel("Frequency (Hz)")
	ylabel("Sₓ(f)")
end


"""
	get_distribution_bound(distribution, p, tail=:center)


Gets the bounds of a normalized distribution where the sum is equal to `p`.


Required Arguments:

- `distribution`: histogram object that contains the bin locations and their
                  weights
- `p`: value rangin from 0 to 1 that specifies the ratio of the distribution to
       include


Optional Arguments:

- `tail`: whether start from the center and obtain bounds above and below the
          center or whether to start from the right or left.
		  `(default = :center)`


Returns:

- `bounds`: 2-element vector that describes the bounds of the distribution which
            include the bins, whose weights sum to equal `p`
"""
function get_distribution_bounds(distribution, p, tail=:center)
	if tail == :center
		center_idx = ceil(Int, length(distribution.edges[1])/2)
		center_idx = findall(x->x==0., distribution.edges[1])[1]
		side_limits = ceil(Int, length(distribution.edges[1])/2)
		distrbution_sum = distribution.weights[center_idx]
		bounds = 0
		for i in 1:side_limits
			if distrbution_sum < p
				distrbution_sum += distribution.weights[center_idx-i]
				distrbution_sum += distribution.weights[center_idx+i]
				bounds += 1
			end
		end
		bounds = [center_idx-bounds, center_idx+bounds]
	elseif tail == :left
		distrbution_sum = 0.
		bounds = [1, 1]
		for i in 1:length(distribution.weights)
			if distrbution_sum < p
				distrbution_sum += distribution.weights[i]
				bounds[2] = i
			end
		end
	elseif tail == :right
		distrbution_sum = 0.
		bounds = [1, 1] .* length(distribution.weights)
		for i in length(distribution.weights):-1:1
			if distrbution_sum < p
				distrbution_sum += distribution.weights[i]
				bounds[1] = i
			end
		end
	else
		error("$(tail) invalid. Must be either `:center`, `:left`, or `:right`.")
	end
	return bounds
end


"""
	snr2cn0(SNR, T)


Convert SNR to C/N₀ for an SNR estimate determined
at an integration time of `T`.


Required Arguments:

- `SNR`: the signal-to-noise ratio in dB
- `T`: the integration time in seconds

Returns:


- `CN0`: carrier-to-noise ratio in dB⋅Hz
"""
function snr2cn0(SNR, T)
	return SNR + 10*log10(1/T)
end


"""
	phase_noise_variance(CN0, T)


Calculate the expected carrier phase measurement variance, `σᵩ²`. 


Required Arguments:

- `CN0`: the carrier-to-noise ratio of the signal in dB⋅Hz
- `T`: the integration time in seconds


Returns:

- `σᵩ²`: the expected carrier phase measurement variance in rad²
"""
function phase_noise_variance(CN0, T)
	CN0_linear = 10^(CN0/10)
	return (1 + 1/(2*T*CN0_linear)) / (2*T*CN0_linear)
end