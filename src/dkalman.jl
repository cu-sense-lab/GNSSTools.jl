""""
    `dare(A, B, Q, R)`


Compute `X`, the solution to the discrete-time algebraic Riccati equation,
defined as A'XA - X - (A'XB)(B'XB + R)^-1(B'XA) + Q = 0, where Q>=0 and R>0


Algorithm taken from:
Laub, "A Schur Method for Solving Algebraic Riccati Equations."
[http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf](http://dspace.mit.edu/bitstream/handle/1721.1/1301/R-0859-05666488.pdf)


Function taken from the official `ControlSystems.jl` repository
[https://github.com/JuliaControl/ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl)
"""
function dare(A, B, Q, R)
    if (!ishermitian(Q) || minimum(eigvals(real(Q))) < 0)
        error("Q must be positive-semidefinite.");
    end
    if (!isposdef(R))
        error("R must be positive definite.");
    end

    n = size(A, 1);

    E = [
        Matrix{Float64}(I, n, n) B/R*B';
        zeros(size(A)) A'
    ];
    F = [
        A zeros(size(A));
        -Q Matrix{Float64}(I, n, n)
    ];

    QZ = schur(F, E);
    QZ = ordschur(QZ, abs.(QZ.alpha./QZ.beta) .< 1);

    return QZ.Z[(n+1):end, 1:n]/QZ.Z[1:n, 1:n];
end


"""
    `dlqr(A, B, Q, R)`, `dlqr(sys, Q, R)`


Calculate the optimal gain matrix `K` for the state-feedback law `u[k] = K*x[k]` that
minimizes the cost function:

J = sum(x'Qx + u'Ru, 0, inf).

For the discrte time model `x[k+1] = Ax[k] + Bu[k]`. See also `lqg`.


Usage example:

```julia
using LinearAlgebra # For identity matrix I
h = 0.1
A = [1 h; 0 1]
B = [0;1]
C = [1 0]
sys = ss(A,B,C,0, h)
Q = I
R = I
L = dlqr(A,B,Q,R) # lqr(sys,Q,R) can also be used
u(x,t) = -L*x # Form control law,
t=0:h:5
x0 = [1,0]
y, t, x, uout = lsim(sys,u,t,x0=x0)
plot(t,x, lab=["Position"  "Velocity"], xlabel="Time [s]")
```


Function taken from the official `ControlSystems.jl` repository
[https://github.com/JuliaControl/ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl)
"""
function dlqr(A, B, Q, R)
    S = dare(A, B, Q, R)
    K = (B'*S*B + R)\(B'S*A)
    return K
end


"""
    `dkalman(A, C, R1, R2)` kalman(sys, R1, R2)`


Calculate the optimal Kalman gain for discrete time systems.


Function taken from the official `ControlSystems.jl` repository
[https://github.com/JuliaControl/ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl)
"""
dkalman(A, C, R1,R2) = Matrix(dlqr(A',C',R1,R2)')


"""
`place(A, B, p)`, `place(sys::StateSpace, p)`


Calculate gain matrix `K` such that the poles of `(A-BK)` in are in `p`. Uses
Ackermann's formula.


Function taken from the official `ControlSystems.jl` repository
[https://github.com/JuliaControl/ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl)
"""
function place(A, B, p)
    n = length(p)
    n != size(A,1) && error("Must define as many poles as states")
    n != size(B,1) && error("A and B must have same number of rows")
    if size(B,2) == 1
        acker(A,B,p)
    else
        error("place only implemented for SISO systems")
    end
end
