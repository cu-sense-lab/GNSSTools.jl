# """
#     calc_f_code_d(f_d, chipping_rate, sig_freq)
# """
# function calc_f_code_d!(du, u, parms, t)
#     sig_freq, chipping_rate = parms
#     n, f_d, f_d_rate = u
#     f_code_d = chipping_rate*(1 + f_d/sig_freq)
#     f_code_dd = chipping_rate*f_d_rate/sig_freq
#     du[1] = f_code_d + f_code_dd*t
#     du[2] = f_code_dd
#     du[3] = 0.
# end

# """
#     doppler2chips(orbs, chipping_rate, chip_init, sig_freq, t, N,
#                   start_t)
# """
# function doppler2chips(doppler_curve, chipping_rate, chip_init, sig_freq, t;
#                        Δt=t[2]-t[1], N=length(doppler_curve))
#     code_chips = Array{Float64}(undef, N)
#     code_chips = zeros(N)
#     code_chips[1] = chip_init
#     for i in 2:N
#         # x = [code_chips[i-1], doppler_curve[i-1], doppler_curve_diff[i-1]]
#         n = code_chips[i-1]
#         f_d = doppler_curve[i-1]
#         f_d_rate = doppler_curve_diff[i-1]
#         f_d_rate = (doppler_curve[i] - doppler_curve[i-1])/Δt
#         f_code_d = chipping_rate*(1 + f_d/sig_freq)
#         f_code_dd = chipping_rate*f_d_rate/sig_freq
#         # prob = ODEProblem(calc_f_code_d!, x, (0, Δt), parms)
#         # sol = solve(prob, ode45())
#         # code_chips[i] = sol[end][1]
#         code_chips[i] = n + f_code_d*Δt + 0.5*f_code_dd*Δt^2
#     end
#     return code_chips
# end


"""
doppler2chips(signal::ReplicaSignals, doppler_curve,
              doppler_t; Δt=doppler_t[2]-doppler_t[1],
              N=length(doppler_curve))


Use the Doppler curve to generate functions that calculate the code chip and
carrier phase as a function of time.

The Doppler curve is integrated over each `Δt` to deterimine the expected code
phase for a given signal code layer. It is then interpolated and placed into a
vector. There are two separate vectors with interpolated functions inside. One
for the I channel and one for the Q channel. These functions are in the order
they are stored in `signal.code_type.[I or Q]_codes`.


Required Arguments:

- `signal::ReplicaSignals`: replica signal struct
- `doppler_curve`: Doppler frequency curve in Hz
- `doppler_t`: time vector in seconds, corresponding to `doppler_curve`


Optional Arguments:

- `Δt`: time step size in seconds
- `N`: number of samples in `doppler_curve`


Returns:

- `code_chip_I_sitp`: vector of interpolated functions that are functions of `t`
                      for the I channel
- `code_chip_Q_sitp`: vector of interpolated functions that are functions of `t`
                      for the I channel
- `ϕs_sitp`: interpolated function that takes `t` and returns the phase in rads
"""
function doppler2chips(signal::ReplicaSignals, doppler_curve,
                       doppler_t; Δt=doppler_t[2]-doppler_t[1],
                       N=length(doppler_curve))
      # Get parameters from `signal` struct
      ϕ_init = signal.phi
      f_if = signal.f_if
      sig_freq = signal.signal_type.sig_freq
      include_I = signal.signal_type.include_I
      include_Q = signal.signal_type.include_Q
      # chipping_rates_I = signal.signal_type.I_codes.chipping_rates
      # chipping_rates_Q = signal.signal_type.Q_codes.chipping_rates
      # chip_init_I = signal.init_code_phases_I
      # chip_init_Q = signal.init_code_phases_Q
      # Generate range of time for use in interpolation
      t_range = range(doppler_t[1], doppler_t[end]; length=length(doppler_t))
      # code_chips_I = zeros(N, length(signal.signal_type.I_codes.code_num))
      # code_chips_Q = zeros(N, length(signal.signal_type.Q_codes.code_num))
      f_code_d_I_sitps = []
      f_code_d_Q_sitps = []
      # I channel
      if include_I
          chipping_rates_I = signal.signal_type.I_codes.chipping_rates
          chip_init_I = signal.init_code_phases_I
          code_chips_I = zeros(N, signal.signal_type.I_codes.code_num)
          for i in 1:signal.signal_type.I_codes.code_num
              chipping_rate = chipping_rates_I[i]
              code_chips_I[1,i] = chip_init_I[i]
              f_code_d = chipping_rate .* (1 .+ doppler_curve./sig_freq)
              push!(f_code_d_I_sitps, CubicSplineInterpolation(t_range, f_code_d))
          end
      end
      # Q channel
      if include_Q
          chipping_rates_Q = signal.signal_type.Q_codes.chipping_rates
          chip_init_Q = signal.init_code_phases_Q
          code_chips_Q = zeros(N, signal.signal_type.Q_codes.code_num)
          for i in 1:signal.signal_type.Q_codes.code_num
              chipping_rate = chipping_rates_Q[i]
              code_chips_Q[1,i] = chip_init_Q[i]
              f_code_d = chipping_rate .* (1 .+ doppler_curve./sig_freq)
              push!(f_code_d_Q_sitps, CubicSplineInterpolation(t_range, f_code_d))
          end
      end
      ϕs = zeros(N)
      ϕs[1] = ϕ_init
      doppler_curve_rad = (f_if .+ doppler_curve).*2π
      doppler_sitp = CubicSplineInterpolation(t_range, doppler_curve_rad)
      for i in 2:N
          # I Channel
          if include_I
              for j in 1:signal.signal_type.I_codes.code_num
                  code_chips_I[i,j] = code_chips_I[i-1,j] + quadgk(f_code_d_I_sitps[j],
                                                                   (i-2)*Δt+doppler_t[1],
                                                                   (i-1)*Δt+doppler_t[1])[1]
              end
          end
          # Q Channel
          if include_Q
              for j in 1:signal.signal_type.Q_codes.code_num
                  code_chips_Q[i,j] = code_chips_Q[i-1,j] + quadgk(f_code_d_Q_sitps[j],
                                                                   (i-2)*Δt+doppler_t[1],
                                                                   (i-1)*Δt+doppler_t[1])[1]
              end
          end
          # Carrier phase
          ϕs[i] = ϕs[i-1] + quadgk(doppler_sitp, (i-2)*Δt+doppler_t[1], (i-1)*Δt+doppler_t[1])[1]
      end
      # Interpolate I channel codes
      if include_I
          code_chip_I_sitp = Array{typeof(doppler_sitp)}(undef, signal.signal_type.I_codes.code_num)
          for i in 1:signal.signal_type.I_codes.code_num
              code_chip_I_sitp[i] = CubicSplineInterpolation(t_range, view(code_chips_I, :, i))
          end
      else
          code_chip_I_sitp = missing
      end
      # Interpolate Q channel codes
      if include_Q
          code_chip_Q_sitp = Array{typeof(doppler_sitp)}(undef, signal.signal_type.Q_codes.code_num)
          for i in 1:signal.signal_type.Q_codes.code_num
              code_chip_Q_sitp[i] = CubicSplineInterpolation(t_range, view(code_chips_Q, :, i))
          end
      else
          code_chip_Q_sitp = missing
      end
      # Interpolate carrier phase
      ϕs_sitp = CubicSplineInterpolation(t_range, ϕs)
      return (code_chip_I_sitp, code_chip_Q_sitp, ϕs_sitp)
end


"""
    get_chips_and_ϕ(signal::ReplicaSignals, doppler_curve, doppler_t)


Interpolate the Doppler frequency curve and return two functions that calculate
the code and carrier phase, respectively, as a function of time.


Required Arguments:

- `signal::ReplicaSignals`: replica signal struct
- `doppler_curve`: Doppler frequency curve in Hz
- `doppler_t`: time vector in seconds, corresponding to `doppler_curve`


Returns:

- `get_code_val`: function of `t` that returns the current code value for both
                  the I and Q channels
- `get_ϕ`: interpolated function that takes `t` and returns the phase in rads
"""
function get_chips_and_ϕ(signal::ReplicaSignals, doppler_curve, doppler_t)
    code_chips_I, code_chips_Q, get_ϕ = doppler2chips(signal, doppler_curve,
                                                      doppler_t)
    get_code_val(t) = calc_code_val(signal, t, code_chips_I, code_chips_Q)
    return (get_code_val, get_ϕ)
end
