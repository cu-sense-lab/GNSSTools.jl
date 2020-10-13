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
doppler2chips(doppler_curve, chipping_rates, sig_freq, t;
              Δt=t[2]-t[1], N=length(doppler_curve),
              chip_init=zeros(length(chipping_rates)),
              ϕ_init=0.)

`M` should be an integer multiple of `N`
"""
function doppler2chips(doppler_curve, chipping_rates, sig_freq, f_if, t;
                       Δt=t[2]-t[1], N=length(doppler_curve),
                       chip_init=zeros(length(chipping_rates)),
                       ϕ_init=0.)
      t_range = range(t[1], t[end]; length=length(t))
      code_chips = zeros(N, length(chipping_rates))
      f_code_d_sitps = []
      for i in 1:length(chipping_rates)
          chipping_rate = chipping_rates[i]
          code_chips[1,i] = chip_init[i]
          f_code_d = chipping_rate .* (1 .+ doppler_curve./sig_freq)
          push!(f_code_d_sitps, CubicSplineInterpolation(t_range, f_code_d))
      end
      ϕs = zeros(N)
      ϕs[1] = ϕ_init
      doppler_curve_rad = (f_if .+ doppler_curve).*2π
      doppler_sitp = CubicSplineInterpolation(t_range, doppler_curve_rad)
      for i in 2:N
          for j in 1:length(chipping_rates)
              code_chips[i,j] = code_chips[i-1,j] + quadgk(f_code_d_sitps[j],
                                                           (i-2)*Δt+t[1],
                                                           (i-1)*Δt+t[1])[1]
          end
          ϕs[i] = ϕs[i-1] + quadgk(doppler_sitp, (i-2)*Δt+t[1], (i-1)*Δt+t[1])[1]
      end
      code_chip_sitp = []
      for i in 1:length(chipping_rates)
          push!(code_chip_sitp, CubicSplineInterpolation(t_range, view(code_chips, :, i)))
      end
      ϕs_sitp = CubicSplineInterpolation(t_range, ϕs)
      return (code_chip_sitp, ϕs_sitp)
end


"""
    get_chips_and_ϕ(signal::L1CASignal, doppler_curve, doppler_t)
"""
function get_chips_and_ϕ(signal::L1CASignal, doppler_curve; doppler_t=signal.t)
    sig_freq = signal.sig_freq
    f_if = signal.f_if
    primary_chipping_rate = signal.chipping_rate
    databit_chipping_rate = l1ca_db_chipping_rate
    chipping_rates = [primary_chipping_rate, databit_chipping_rate]
    primary_code_init = signal.l1ca_init_code_phase
    databit_code_init = signal.db_init_code_phase
    chip_inits = [primary_code_init, databit_code_init]
    ϕ_init = signal.ϕ
    code_chips, get_ϕ = doppler2chips(doppler_curve, chipping_rates, sig_freq,
                                      f_if, doppler_t; chip_init=chip_inits,
                                      ϕ_init=ϕ_init)
    get_code_val(t) = calc_code_val(signal, t, code_chips)
    return (get_code_val, get_ϕ)
end


"""
    get_chips_and_ϕ(signal::L5QSignal, doppler_curve, doppler_t)
"""
function get_chips_and_ϕ(signal::L5QSignal, doppler_curve; doppler_t=signal.t)
    sig_freq = signal.sig_freq
    f_if = signal.f_if
    primary_chipping_rate = signal.chipping_rate
    seconday_chipping_rate = nh_chipping_rate
    chipping_rates = [primary_chipping_rate, seconday_chipping_rate]
    primary_code_init = signal.l5q_init_code_phase
    secondary_code_init = signal.nh_init_code_phase
    chip_inits = [primary_code_init, secondary_code_init]
    ϕ_init = signal.ϕ
    code_chips, get_ϕ = doppler2chips(doppler_curve, chipping_rates, sig_freq,
                                      f_if, doppler_t; chip_init=chip_inits,
                                      ϕ_init=ϕ_init)
    get_code_val(t) = calc_code_val(signal, t, code_chips)
    return (get_code_val, get_ϕ)
end


"""
    get_chips_and_ϕ(signal::L5ISignal, doppler_curve, doppler_t)
"""
function get_chips_and_ϕ(signal::L5ISignal, doppler_curve; doppler_t=signal.t)
    sig_freq = signal.sig_freq
    f_if = signal.f_if
    primary_chipping_rate = signal.chipping_rate
    seconday_chipping_rate = nh_chipping_rate
    databit_chipping_rate = L5_db_chipping_rate
    chipping_rates = [primary_chipping_rate, seconday_chipping_rate, databit_chipping_rate]
    primary_code_init = signal.l5i_init_code_phase
    secondary_code_init = signal.nh_init_code_phase
    databit_code_init = signal.db_init_code_phase
    chip_inits = [primary_code_init, secondary_code_init, databit_code_init]
    ϕ_init = signal.ϕ
    code_chips, get_ϕ = doppler2chips(doppler_curve, chipping_rates, sig_freq,
                                      f_if, doppler_t; chip_init=chip_inits,
                                      ϕ_init=ϕ_init)
    get_code_val(t) = calc_code_val(signal, t, code_chips)
    return (get_code_val, get_ϕ)
end
