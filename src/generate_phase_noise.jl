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
    generate_phase_noise(t_length, f_s, h₋₂=2*4e-22, h₋₁=2*4e-21, h₀=3*9e-22,
                         h₁=2*5e-23, h₂=2*0e-26)


Generate phase noise using the h-paramters `hₐ` where `a` is [-2, -1, 0, 1, 2].
"""
function generate_phase_noise2(t_length, f_s, v_0=f_s, h₋₂=2*4e-22, h₋₁=2*4e-21, h₀=3*9e-22,
                              h₁=2*5e-23, h₂=2*0e-26)
    N = floor(Int, t_length*f_s)
    N_over_2 = floor(Int, N/2)
    # Square root of the spectral density
    Δf = 1/t_length
    S_y_sqrt = zeros(N)
    for i in 1:(N_over_2+1)
        f = (i-1)*Δf
        if i == 1
            val = 0
        else
            val = v_0^2 * (h₋₂*f^(-2) + h₋₁*f^(-1) + h₀*f^0 + h₁*f^1 + h₂*f^2)
            # val = v_0^2 * (h₋₂*f^(-4) + h₋₁*f^(-3) + h₀*f^(-2) + h₁*f^(-1) + h₂*f^0)
        end
        if (i > 1) && (i < (N_over_2+1))
            val = val/2
        end
        val = sqrt(val*f_s*N)
        S_y_sqrt[i] = val
    end
    # Copy the positive frequency powers to the negative frequencies, except
    # the DC and Nyquist frequency
    S_y_sqrt[N_over_2+2:end] = reverse(S_y_sqrt[2:N_over_2])
    # return S_y_sqrt
    # return ifft(S_y_sqrt.*randn(Complex{Float64}, N))
    return ifft(S_y_sqrt.*fft(randn(N)))
end


function generate_phase_noise(t_length, f_s; hₐ=1, α=-1)
    N = floor(Int, t_length*f_s)
    T₀ = 1/f_s
    Q_d = hₐ/(2 * (2π)^α * T₀^(α-1))
    w = sqrt(Q_d) .* randn(Float64, N)
    h = Array{Float64}(undef, N)
    β = -α/2
    h[1] = 1
    for i in 2:N
        h[i] = h[i-1]*(β + (i-2))/(i-1)
    end
    W = fft(w)
    H = fft(h)
    result = ifft(H.*W)./N
    return (result, h, w, Q_d)
end


function phase_noise(t_length, f_s; h₋₂=2*4e-22, h₋₁=2*4e-21, h₀=3*9e-22)
    N = floor(Int, t_length*f_s)
    τ = 1/f_s
    Q_d_2 = sqrt(h₋₂ * (2π)^2 * τ / 6)
    Q_d_1 = sqrt(h₋₁ * 2 * log(2))
    Q_d_0 = sqrt(h₀ / (2*τ))
    w_2 = Q_d_2 .* randn(Float64, N)
    w_1 = Q_d_1 .* randn(Float64, N)
    w_0 = Q_d_0 .* randn(Float64, N)
    h_2 = zeros(N)
    h_1 = zeros(N)
    h_0 = zeros(N)
    h_2[1] = 1
    h_1[1] = 1
    h_0[1] = 1
    for i in 2:N
        h_2[i] = (i - 2 - (-1.99/2)) * h_2[i-1] / (i-1)
        h_1[i] = (i - 2 - (-1/2)) * h_1[i-1] / (i-1)
        h_0[i] = (i - 2 - (-0/2)) * h_0[i-1] / (i-1)
    end
    H_2 = fft(h_2)
    H_1 = fft(h_1)
    H_0 = fft(h_0)
    W_2 = fft(w_2)
    W_1 = fft(w_1)
    W_0 = fft(w_0)
    w_2 = ifft(W_2 .* H_2)
    w_1 = ifft(W_1 .* H_1)
    w_0 = ifft(W_0 .* H_0)
    w = w_2.+w_1.+w_0
    return (w_2, w_1, w_0, w)
end


function test(t_length, f_s; hₐ=1, α=-1, Nscale=2, spectrum=false, log_X=true, figsize=missing)
    results, h, w, Q_d = generate_phase_noise(t_length, f_s; hₐ=hₐ, α=α); 
    if ismissing(figsize)
        figure();
    else
        figure(figsize=figsize);
    end 
    subplot(1,2,1);
    if spectrum
        plot_spectrum(results, f_s; log_freq=true, new_fig=false, log_X=log_X)
    else
        psd(results, Fs=f_s, NFFT=Int(length(results)/Nscale)); 
    end
    xscale("log"); 
    subplot(1,2,2); 
    plot(Array(range(0, t_length, length=length(results))), results);
    xlabel("Time (s)")
    ylabel("x(t)")
    return (results, h, w, Q_d)
end
