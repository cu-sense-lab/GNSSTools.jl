"""
    calc_code_val(signal::L5QSignal, t)

Calculates the value of the L5Q code with
parameters defined by `signal` at a time
specified by `t` in seconds. Returns an
Int64 value that is either -1 or 1.
"""
function calc_code_val(signal::L5QSignal, t)
    # Get L5Q code value at t
    l5q = l5q_codes[signal.prn][calccodeidx(signal.l5q_init_code_phase,
                                            signal.f_l5q_d, signal.f_l5q_dd,
                                            t, L5_code_length)]
    # Get Neuman code sequence value at t
    nh = nh20[calccodeidx(signal.nh_init_code_phase,
                          signal.f_nh_d, signal.f_nh_dd,
                          t, nh_code_length)]
    return xor(l5q, nh)*2-1
end
