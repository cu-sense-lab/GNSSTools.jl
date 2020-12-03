"""
    generate_phase_noise(N; scale=1/10)
"""
function generate_phase_noise(N; scale=1/10)
    noise = randn(N)
    noise_temp = randn(N)
    n = 2
    sampled_noise_T = Array{Int}(undef, 0)
    while n <= N
        push!(sampled_noise_T, n)
        n = 2*n
    end
    M = length(sampled_noise_T)
    idx = 1
    T = 2
    for i in 1:N
        if idx > N
            idx = 1
            T = 2*T
        end
        if (idx+T-1) > N
            view(noise, idx:N) .+= noise_temp[i]
        else
            view(noise, idx:idx+T-1) .+= noise_temp[i]
        end
        idx += T
    end
    scale_factor = scale/M
    view(noise, 1:N) .*= scale_factor
    return noise
end


"""
    generate_phase_noise!(noise, N=length(noise); scale=1/10)
"""
function generate_phase_noise!(noise, N=length(noise)::Int; scale=1/10)
    noise_temp = randn(N)
    n = 2
    sampled_noise_T = Array{Int}(undef, 0)
    while n <= N
        push!(sampled_noise_T, n)
        n = 2*n
    end
    M = length(sampled_noise_T)
    idx = 1
    T = 2
    for i in 1:N
        if idx > N
            idx = 1
            T = 2*T
        end
        if (idx+T-1) > N
            view(noise, idx:N) .+= noise_temp[i]
        else
            view(noise, idx:idx+T-1) .+= noise_temp[i]
        end
        idx += T
    end
    scale_factor = scale/M
    view(noise, 1:N) .*= scale_factor
    return noise
end


"""
    generate_phase_noise!(noise, temp_noise; scale=1/10)
"""
function generate_phase_noise!(noise, noise_temp::Vector; scale=1/10)
    n = 2
    sampled_noise_T = Array{Int}(undef, 0)
    N = length(noise)
    while n <= N
        push!(sampled_noise_T, n)
        n = 2*n
    end
    M = length(sampled_noise_T)
    idx = 1
    T = 2
    for i in 1:N
        if idx > N
            idx = 1
            T = 2*T
        end
        if (idx+T-1) > N
            view(noise, idx:N) .+= noise_temp[i]
        else
            view(noise, idx:idx+T-1) .+= noise_temp[i]
        end
        idx += T
    end
    scale_factor = scale/M
    view(noise, 1:N) .*= scale_factor
    return noise
end


"""
    plot_spectrum(x)
"""
function plot_spectrum(x, f_s, flog=false)
    N = length(x)
    t_length = N/f_s
    freqs = f = Array(0:1/t_length:f_s/2)
    X = 20*log10.(abs2.(fft(x)))[1:length(freqs)]
    fig = figure()
    ax = fig.add_subplot(1, 1, 1)
    if flog
        ax.set_xscale("log")
    end
    plot(freqs, X)
end
