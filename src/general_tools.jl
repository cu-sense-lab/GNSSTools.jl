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
    gnsstypes

Dictionary containing the qyuivalent strings for each
type used in `GNSSTools`.
"""
const gnsstypes = Dict(Val{:l5q}() => "l5q",
                       Val{:l5i}() => "l5i",
                       Val{:l1ca}() => "l1ca",
                       Val{:fft}() => "fft",
                       Val{:carrier}() => "carrier",
                       Val{:sc8}() => "sc8",
                       Val{:sc4}() => "sc4")


"""
    calcinitcodephase(code_length, f_code_d, f_code_dd,
                      f_s, code_start_idx)

Calculates the initial code phase of a given code
where f_d and fd_rate are the Doppler affected
code frequency and code frequency rate, respectively.
"""
function calcinitcodephase(code_length, f_code_d, f_code_dd,
                           f_s, code_start_idx)
    t₀ = (code_start_idx-1)/f_s
    init_phase = -f_code_d*t₀ - 0.5*f_code_dd*t₀^2
    return (init_phase%code_length + code_length)%code_length
end


"""
    calccodeidx(init_chip, f_code_d, f_code_dd, t, code_length)

Calculates the index in the codes for a given t.
"""
function calccodeidx(init_chip, f_code_d, f_code_dd,
                     t, code_length)
    return Int(floor(init_chip+f_code_d*t+0.5*f_code_dd*t^2)%code_length)+1
end


"""
    calctvector(N, f_s)

Generates a `N` long time vector
with time spacing `Δt` or `1/f_s`.
"""
function calctvector(N, f_s)
    # Generate time vector
    t = Array{Float64}(undef, N)
    @threads for i in 1:N
        @inbounds t[i] = (i-1)/f_s
    end
    return t
end


"""
    meshgrid(x, y)

Generate a meshgrid the way Python would in Numpy.
"""
function meshgrid(x, y)
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

Find a number sequence in a file name. `digit_nums` can
be a number or array.`separators` can also be either a
number or array containing the characters such as `_` and `.`
that may be between the numbers. Returns the string containing
the sequence if found. Returns `missing` if not found.
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

Find a sequency of 8 digits with `_` separating it from a sequence
of 6 digits. Return the timestamp tuple.
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
    get_signal_type(file_name)

Find the signal type of the data based off its file name.
Only checks if signal type is L1 or L5 signal. L2 is not
supported.
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
"""
function make_subplot(fig, row, col, i; projection3d=false, aspect="auto")
	if projection3d
		mplot3d = PyPlot.PyObject(PyPlot.axes3D)  # must be called in local scope
		                                          # in order to make 3D subplots
		ax = fig.add_subplot(row, col, i, projection="3d")
	else
		ax = fig.add_subplot(row, col, i, aspect=aspect)
	end
	return (fig, ax)
end


"""
	calc_doppler_code_rate(f_code, f_carrier, f_d)
"""
function calc_doppler_code_rate(f_code, f_carrier, f_d, fd_rate)
	f_code_d = f_code*(1. + f_d/f_carrier)
	f_code_dd = f_code*fd_rate/f_carrier
	return (f_code_d, f_code_dd)
end
