"""
    voss!(noise1::Vector, noise2::Vector, scale, N=length(noise1))


Generates pink noise using two random Gaussian noise vectors and the
Voss-McCartney algorithm. It is described in detaile at
[https://www.firstpr.com.au/dsp/pink-noise/#Voss-McCartney](https://www.firstpr.com.au/dsp/pink-noise/#Voss-McCartney).


Required Arguments:

- `noise1::Vector`: length `N` Gaussian noise vector
- `noise2::Vector`: length `N` Gaussian noise vector
- `scale`: the amplitude of the noise in rad


Optional Arguments:

- `N`: length of the noise vectors `(default = length(noise1))`


Modifies and Returns:

- `noise1`: length `N` pink noise vector
"""
function voss!(noise1::Vector, noise2::Vector, scale, N=length(noise1))
    # `n` is the starting period of the each noise vector, the number of samples
    # before the noise value changes for a given layer
    n = 2
    # For a given `N`, determine all the periods of the noise that are a
    # multiple of `n`. Each noise layer has a period that is twice the period
    # from the previous layer, until the period is equal to the longest, `N`.
    # `sampled_noise_T` stores all the periods that are that used in
    # constructing pink noise.
    sampled_noise_T = Array{Int}(undef, 0)
    while n <= N
        push!(sampled_noise_T, n)
        n = 2*n
    end
    # `noise1` vector is noise with a period of 1 sample. The second vector,
    # `noise2`, contains all noise layers. The first half of `noise2` is the
    # second layer of the pink noise, where the period is 2 samples. The rest
    # of `noise2` is used all other periods that are multiples of 2 until the
    # period is equal to `N`.
    M = length(sampled_noise_T)
    idx = 1
    T = 2
    for i in 1:N
        if idx > N
            idx = 1
            T = 2*T
        end
        if (idx+T-1) > N
            view(noise1, idx:N) .+= noise2[i]
        else
            view(noise1, idx:idx+T-1) .+= noise2[i]
        end
        idx += T
    end
    # Scale the amplitude of `noise1` by the number of layers and `scale`.
    scale_factor = scale/M
    view(noise1, 1:N) .*= scale_factor
    return noise1
end


"""
    generate_phase_noise(N::Int; scale=1/10)


Generates a pink noise vector of length `N` with noise σ of `scale` rad. Two
Gaussian noise vectors are generated and passed to `voss!`. Pink noise is
stored in `noise`.


Required Arguments:

- `N::Int`: the length of the noise vector


Optional Arguments:

- `scale`: the amplitude of the noise in rad `(default = 1/10)`


Returns:

- `noise::Vector`: contains pink noise with amplitude `scale` radians
"""
function generate_phase_noise(N::Int; scale=1/10)
    noise = randn(N)
    noise_temp = randn(N)
    return voss!(noise, noise_temp, scale, N)
end


"""
    generate_phase_noise!(noise::Vector, N::Int=length(noise); scale=1/10)


Generates pink noise and stores it in `noise`. A second noise vector of length
`N` is generated and both noise vectors are passed to `voss!`.


Required Arguments:

- `noise::Vector`: length `N` Gaussian noise vector


Optional Arguments:

- `N::Int`: the length of the noise vector `(default = length(noise))`
- `scale`: the amplitude of the noise in rad `(default = 1/10)`


Returns:

- `noise::Vector`: contains pink noise with amplitude `scale` radians
"""
function generate_phase_noise!(noise::Vector, N::Int=length(noise); scale=1/10)
    noise_temp = randn(N)
    return voss!(noise, noise_temp, scale, N)
end


"""
    generate_phase_noise!(noise::Vector, noise_temp::Vector; scale=1/10)


Generates pink noise using two Gaussian noise generated vectors which are passed
to `voss!`. Pink noise is generated and stored in the `noise` vector.


Required Arguments:

- `noise::Vector`: length `N` Gaussian noise vector
- `noise_temp::Vector`: length `N` Gaussian noise vector


Optional Arguments:

- `scale`: the amplitude of the noise in rad `(default = 1/10)`


Returns:

- `noise::Vector`: contains pink noise with amplitude `scale` radians
"""
function generate_phase_noise!(noise::Vector, noise_temp::Vector; scale=1/10)
    return voss!(noise, noise_temp, scale, length(noise))
end


"""
    generate_phase_noise(t_length::Number, f_s::Number, f_carrier::Number, 
                         h_parms::Vector)


Generate phase noise using the h-paramters, stored in the vector `h_parms` where
the order is `[h₋₂, h₋₁, h₀, h₁, h₂]`.


Required Arguments:

- `t_length`: length of of the phase noise array to generate in seconds
- `f_s`: sampling rate of the data in Hz
- `f_carrier`: receiver center frequency in Hz (i.e. L1 C/A carrier frequency)
- `h_parms`: vector containing the oscillator h parameters in the from
             of `[h₋₂, h₋₁, h₀, h₁, h₂]`


Returns:

- `phase_noise`: vector of generated oscillator phase noise in radians
"""
function generate_phase_noise(t_length::Number, f_s::Number, f_carrier::Number, 
                              h_parms::Vector)
    N = floor(Int, t_length*f_s)
    phase_noise = Array{Complex{Float64}}(undef, N)
    generate_phase_noise!(phase_noise, t_length,  f_carrier, h_parms)
    return phase_noise
end


"""
    generate_phase_noise!(phase_noise::Vector t_length::Number,
                          f_carrier::Number, h_parms::Vector)


Generate phase noise using the h-paramters, stored in the vector `h_parms` where
the order is `[h₋₂, h₋₁, h₀, h₁, h₂]`. Performs in-place operation on 
pre-existing `phase_noise` vector.

Generates phase noise using methods from
[https://www.mathworks.com/matlabcentral/fileexchange/47342-timeseriesfrompsd-sxx-fs-t-plot_on](https://www.mathworks.com/matlabcentral/fileexchange/47342-timeseriesfrompsd-sxx-fs-t-plot_on)


Required Arguments:

Required Arguments:

- `phase_noise`: vector to store generated oscillator phase noise 
- `t_length`: length of of the phase noise array to generate in seconds
- `f_carrier`: receiver center frequency in Hz (i.e. L1 C/A carrier frequency)
- `h_parms`: vector containing the oscillator h parameters in the from
             of `[h₋₂, h₋₁, h₀, h₁, h₂]`


Returns:

- `phase_noise`: vector of generated oscillator phase noise in radians
"""
function generate_phase_noise!(phase_noise::Vector, t_length::Number, 
                               f_carrier::Number, h_parms::Vector)
    h₋₂, h₋₁, h₀, h₁, h₂ = h_parms
    N = length(phase_noise)
    N_over_2 = floor(Int, N/2)
    rand!(view(phase_noise, 1:N_over_2+1), Float64)
    # Square root of the spectral density
    Δf = 1/t_length
    # f = Δf
    # val = f_carrier^2 * (h₋₂*f^(-4) + h₋₁*f^(-3) + h₀*f^(-2) + h₁*f^(-1) + h₂*f^0)
    val = 0
    phase_noise[1] = val * cis(2π * phase_noise[1])
    for i in 2:(N_over_2+1)
        f = (i-1)*Δf
        val = f_carrier^2 * (h₋₂*f^(-4) + h₋₁*f^(-3) + h₀*f^(-2) + h₁*f^(-1) + h₂*f^0)
        if (i > 1) && (i < (N_over_2+1))
            val = val/2
        end
        val = sqrt(val) * cis(2π * phase_noise[i])
        phase_noise[i] = val
    end
    # Copy the positive frequency powers to the negative frequencies, except
    # the DC and Nyquist frequency
    phase_noise[N_over_2+2:end] .= reverse(conj.(phase_noise[2:N_over_2]))
    ifft!(phase_noise)
    for i in 1:N
        phase_noise[i] = phase_noise[i] * N * Δf
    end
    return phase_noise
end
