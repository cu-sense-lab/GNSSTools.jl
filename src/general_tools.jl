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
                       Val{:fft}() => "fft",
                       Val{:carrier}() => "carrier",
                       Val{:sc8}() => "sc8",
                       Val{:sc4}() => "sc4")
