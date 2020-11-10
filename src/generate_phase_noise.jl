"""
    generate_phase_noise(m, N)
"""
function generate_phase_noise(N; scale=1/10)
    noise = randn(N)
    N₀ = 2
    sampled_noise_T = Array{Int}(undef, 0)
    n = N₀
    while n <= N
        push!(sampled_noise_T, n)
        n = 2*n
    end
    M = length(sampled_noise_T)
    sampled_idx = ones(M)
    sampled_noise = randn(M)
    scale_factor = scale/M
    @threads for i in 1:N
        for j in 1:M
            if i%sampled_noise_T[j] == 0
                sampled_noise[j] = randn()
            end
            noise[i] += sampled_noise[j]
        end
        noise[i] = noise[i]*scale_factor
    end
    return noise
end
# function generate_phase_noise(N; scale=1/100)
#     noise = randn(N)
#     noise_temp = Array{Float64}(undef, ceil(Int, N/2))
#     repeat_temp = Array{Float64}(undef, ceil(Int, 2*N))
#     n = 2
#     sampled_noise_T = Array{Int}(undef, 0)
#     while n <= N
#         push!(sampled_noise_T, n)
#         n = 2*n
#     end
#     M = length(sampled_noise_T)
#     @threads for i in 1:M
#         n = sampled_noise_T[i]
#         sample_num = ceil(Int, N/n)
#         randn!(view(noise_temp, 1:sample_num))
#         # println([n, sample_num])
#         view(repeat_temp, 1:n*sample_num) .= repeat(view(noise_temp, 1:sample_num), inner=n)
#         noise .+= view(repeat_temp, 1:N)
#     end
#     scale_factor = scale/M
#     return noise.*scale_factor
# end


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
