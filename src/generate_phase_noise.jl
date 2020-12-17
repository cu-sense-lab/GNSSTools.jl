"""
    voss!(noise1::Vector, noise2::Vector, scale, N=length(noise1))


#


Required Arguments:

- `noise1::Vector`:
- `noise2::Vector`:
- `scale`


Optional Arguments:

- `N`:


Modifies and Returns:

- `noise1`:
"""
function voss!(noise1::Vector, noise2::Vector, scale, N=length(noise1))
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
            view(noise1, idx:N) .+= noise2[i]
        else
            view(noise1, idx:idx+T-1) .+= noise2[i]
        end
        idx += T
    end
    scale_factor = scale/M
    view(noise1, 1:N) .*= scale_factor
    return noise1
end


"""
    generate_phase_noise(N::Int; scale=1/10)


#


Required Arguments:

- `N::Int`:


Optional Arguments:

- `scale`:


Returns:

-
"""
function generate_phase_noise(N::Int; scale=1/10)
    noise = randn(N)
    noise_temp = randn(N)
    return voss!(noise, noise_temp, scale, N)
end


"""
    generate_phase_noise!(noise::Vector, N::Int=length(noise); scale=1/10)


#


Required Arguments:

- `noise::Vector`:


Optional Arguments:

- `N::Int`:
- `scale`:


Returns:

-
"""
function generate_phase_noise!(noise::Vector, N::Int=length(noise); scale=1/10)
    noise_temp = randn(N)
    return voss!(noise, noise_temp, scale, N)
end


"""
    generate_phase_noise!(noise::Vector, noise_temp::Vector; scale=1/10)


#


Required Arguments:

- `noise::Vector`:
- `noise_temp::Vector`:


Optional Arguments:

- `scale`:


Returns:

-
"""
function generate_phase_noise!(noise::Vector, noise_temp::Vector; scale=1/10)
    return voss!(noise, noise_temp, scale, length(noise))
end
