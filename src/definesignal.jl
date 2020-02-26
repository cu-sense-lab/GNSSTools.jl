##################################################################
######################### L1 C/A Signal ##########################
##################################################################
"""
    definesignal(type::Val{:l1ca}, prn, f_s, t_length;
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., ϕ=0., nADC=4, B=2.046e6,
                 include_carrier=true, include_adc=true,
                 include_noise=true, code_start_idx=1)

Define properties of locally generated L1 C/A signal
based off its type, PRN, etc.
"""
function definesignal(type::Val{:l1ca}, f_s, t_length; prn=1,
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=4, B=2.046e6,
                      include_carrier=true, include_adc=true,
                      include_noise=true, code_start_idx=1,
                      include_databits=true)
    sample_num = Int(f_s * t_length)
    # Generate time vector
    t = Array{Float64}(undef, sample_num)
    @threads for i in 1:sample_num
        @inbounds t[i] = (i-1)/f_s  # s
    end
    # Calculate number of data bits (use `ceil` to ensure that databits do not repeat)
    db_code_length = Int64(ceil(l1ca_db_chipping_rate*t_length))
    # Calculate code chipping rates with Doppler applied
    # L1 C/A
    f_l1ca_d = l1ca_chipping_rate*(1. + f_d/L1_freq)
    f_l1ca_dd = l1ca_chipping_rate*fd_rate/L1_freq
    # Data bit sequence
    f_db_d = l1ca_db_chipping_rate*(1. + f_d/L1_freq)
    f_db_dd = l1ca_db_chipping_rate*fd_rate/L1_freq
    # Calculate the L1 C/A and data bit code phase offsets
    l1ca_init_code_phase = calcinitcodephase(l1ca_code_length,
                                             f_l1ca_d, f_l1ca_dd,
                                             f_s, code_start_idx)
    db_init_code_phase = calcinitcodephase(db_code_length,
                                           f_db_d, f_db_dd,
                                           f_s, code_start_idx)
    # Allocate space for signal
    data = Array{Complex{Float64}}(undef, sample_num)
    isreplica = false
    noexp = false
    # Generate random databit vector
    databits = rand(0:1, db_code_length)
    return L1CASignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                      Tsys, CN0, ϕ, nADC, B, code_start_idx,
                      l1ca_init_code_phase, db_init_code_phase, t,
                      data, include_carrier, include_adc,
                      include_noise, include_databits,
                      f_l1ca_d, f_l1ca_dd,
                      f_db_d, f_db_dd, sample_num,
                      isreplica, noexp, l1ca_chipping_rate,
                      L1_freq, l1ca_code_length, db_code_length, databits)
end


"""
    definesignal!(signal::L1CASignal;
                  prn=signal.prn, f_d=signal.f_d,
                  f_if=signal.f_if, fd_rate=signal.fd_rate,
                  Tsys=signal.Tsys, CN0=signal.CN0,
                  ϕ=signal.ϕ, nADC=signal.nADC,
                  B=signal.B,
                  include_carrier=signal.include_carrier,
                  include_adc=signal.include_adc,
                  include_noise=signal.include_noise,
                  code_start_idx=signal.code_start_idx,
                  isreplica=signal.isreplica,
                  include_databits=signal.include_databits)

Redefine properties of locally generated L5Q signal
based off its type, PRN, etc. 

Sampling rate (f_s) and signal length (t_length)
are the only parameters that cannot be redefined.

"""
function definesignal!(signal::L1CASignal;
                       prn=signal.prn, f_d=signal.f_d,
                       f_if=signal.f_if, fd_rate=signal.fd_rate,
                       Tsys=signal.Tsys, CN0=signal.CN0,
                       ϕ=signal.ϕ, nADC=signal.nADC,
                       B=signal.B,
                       include_carrier=signal.include_carrier,
                       include_adc=signal.include_adc,
                       include_noise=signal.include_noise,
                       code_start_idx=signal.code_start_idx,
                       isreplica=signal.isreplica,
                       noexp=signal.noexp,
                       include_databits=signal.include_databits)
    ## Calculate code chipping rates with Doppler applied
    # L1 C/A
    f_l1ca_d = l1ca_chipping_rate*(1. + f_d/L1_freq)
    f_l1ca_dd = l1ca_chipping_rate*fd_rate/L1_freq
    # Data bit sequence
    f_db_d = l1ca_db_chipping_rate*(1. + f_d/L1_freq)
    f_db_dd = l1ca_db_chipping_rate*fd_rate/L1_freq
    # Calculate the L1 C/A and data bit code phase offsets
    l1ca_init_code_phase = calcinitcodephase(l1ca_code_length,
                                            f_l1ca_d, f_l1ca_dd,
                                            signal.f_s,
                                            code_start_idx)
    db_init_code_phase = calcinitcodephase(signal.db_code_length,
                                           f_db_d, f_db_dd,
                                           signal.f_s,
                                           code_start_idx)
    # Store udated variables to `L1CASignal` struct
    signal.prn = prn
    signal.f_d = f_d
    signal.f_if = f_if
    signal.fd_rate = fd_rate
    signal.Tsys = Tsys
    signal.CN0 = CN0
    signal.ϕ = ϕ
    signal.nADC = nADC
    signal.B = B
    signal.include_carrier = include_carrier
    signal.include_adc = include_adc
    signal.include_noise = include_noise
    signal.f_l1ca_d = f_l1ca_d
    signal.f_l1ca_dd = f_l1ca_dd
    signal.f_db_d = f_db_d
    signal.f_db_dd = f_db_dd
    signal.l1ca_init_code_phase = l1ca_init_code_phase
    signal.db_init_code_phase = db_init_code_phase
    signal.code_start_idx = code_start_idx
    signal.isreplica = isreplica
    signal.noexp = noexp
    signal.include_databits = include_databits
    return signal
end


##################################################################
########################### L5Q Signal ###########################
##################################################################
"""
    definesignal(type::Val{:l5q}, prn, f_s, t_length;
                 f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                 CN0=45., ϕ=0., nADC=4, B=2.046e7,
                 include_carrier=true, include_adc=true,
                 include_noise=true, code_start_idx=1)

Define properties of locally generated L5Q signal
based off its type, PRN, etc.
"""
function definesignal(type::Val{:l5q}, f_s, t_length; prn=1,
                      f_if=0., f_d=0., fd_rate=0., Tsys=535.,
                      CN0=45., ϕ=0., nADC=4, B=2.046e7,
                      include_carrier=true, include_adc=true,
                      include_noise=true, code_start_idx=1)
    sample_num = Int(f_s * t_length)
    # Generate time vector
    t = Array{Float64}(undef, sample_num)
    @threads for i in 1:sample_num
        @inbounds t[i] = (i-1)/f_s  # s
    end
    # Calculate code chipping rates with Doppler applied
    # L5Q
    f_l5q_d = L5_chipping_rate*(1. + f_d/L5_freq)
    f_l5q_dd = L5_chipping_rate*fd_rate/L5_freq
    # Neuman sequence
    f_nh_d = nh_chipping_rate*(1. + f_d/L5_freq)
    f_nh_dd = nh_chipping_rate*fd_rate/L5_freq
    # Calculate the L5Q and nh code phase offsets
    l5q_init_code_phase = calcinitcodephase(L5_code_length,
                                            f_l5q_d, f_l5q_dd,
                                            f_s, code_start_idx)
    nh_init_code_phase = calcinitcodephase(nh20_code_length,
                                           f_nh_d, f_nh_dd,
                                           f_s, code_start_idx)
    # Allocate space for signal
    data = Array{Complex{Float64}}(undef, sample_num)
    isreplica = false
    noexp = false
    return L5QSignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
                     Tsys, CN0, ϕ, nADC, B, code_start_idx,
                     l5q_init_code_phase, nh_init_code_phase, t,
                     data, include_carrier, include_adc,
                     include_noise,
                     f_l5q_d, f_l5q_dd,
                     f_nh_d, f_nh_dd, sample_num,
                     isreplica, noexp, L5_chipping_rate,
                     L5_freq, L5_code_length)
end


"""
    definesignal!(signal::L5QSignal;
                  prn=signal.prn, f_d=signal.f_d,
                  f_if=signal.f_if, fd_rate=signal.fd_rate,
                  Tsys=signal.Tsys, CN0=signal.CN0,
                  ϕ=signal.ϕ, nADC=signal.nADC,
                  B=signal.B,
                  include_carrier=signal.include_carrier,
                  include_adc=signal.include_adc,
                  include_noise=signal.include_noise,
                  code_start_idx=signal.code_start_idx,
                  isreplica=signal.isreplica)

Redefine properties of locally generated L5Q signal
based off its type, PRN, etc. 

Sampling rate (f_s) and signal length (t_length)
are the only parameters that cannot be redefined.

"""
function definesignal!(signal::L5QSignal;
                       prn=signal.prn, f_d=signal.f_d,
                       f_if=signal.f_if, fd_rate=signal.fd_rate,
                       Tsys=signal.Tsys, CN0=signal.CN0,
                       ϕ=signal.ϕ, nADC=signal.nADC,
                       B=signal.B,
                       include_carrier=signal.include_carrier,
                       include_adc=signal.include_adc,
                       include_noise=signal.include_noise,
                       code_start_idx=signal.code_start_idx,
                       isreplica=signal.isreplica,
                       noexp=signal.noexp)
    ## Calculate code chipping rates with Doppler applied
    # L5Q
    f_l5q_d = L5_chipping_rate*(1. + f_d/L5_freq)
    f_l5q_dd = L5_chipping_rate*fd_rate/L5_freq
    # Neuman sequence
    f_nh_d = nh_chipping_rate*(1. + f_d/L5_freq)
    f_nh_dd = nh_chipping_rate*fd_rate/L5_freq
    # Calculate the L5Q and nh code phase offsets
    l5q_init_code_phase = calcinitcodephase(L5_code_length,
                                            f_l5q_d, f_l5q_dd,
                                            signal.f_s,
                                            code_start_idx)
    nh_init_code_phase = calcinitcodephase(nh20_code_length,
                                           f_nh_d, f_nh_dd,
                                           signal.f_s,
                                           code_start_idx)
    # Store udated variables to "L5QSignal" struct
    signal.prn = prn
    signal.f_d = f_d
    signal.f_if = f_if
    signal.fd_rate = fd_rate
    signal.Tsys = Tsys
    signal.CN0 = CN0
    signal.ϕ = ϕ
    signal.nADC = nADC
    signal.B = B
    signal.include_carrier = include_carrier
    signal.include_adc = include_adc
    signal.include_noise = include_noise
    signal.f_l5q_d = f_l5q_d
    signal.f_l5q_dd = f_l5q_dd
    signal.f_nh_d = f_nh_d
    signal.f_nh_dd = f_nh_dd
    signal.l5q_init_code_phase = l5q_init_code_phase
    signal.nh_init_code_phase = nh_init_code_phase
    signal.code_start_idx = code_start_idx
    signal.isreplica = isreplica
    signal.noexp = noexp
    return signal
end


# ##################################################################
# ########################### L5I Signal ###########################
# ##################################################################
# """
#     definesignal(type::Val{:l5i}, prn, f_s, t_length;
#                  f_if=0., f_d=0., fd_rate=0., Tsys=535.,
#                  CN0=45., ϕ=0., nADC=4, B=2.046e7,
#                  include_carrier=true, include_adc=true,
#                  include_noise=true, code_start_idx=1,
#                  include_databits=true)

# Define properties of locally generated L5Q signal
# based off its type, PRN, etc.
# """
# function definesignal(type::Val{:l5i}, f_s, t_length; prn=1,
#                       f_if=0., f_d=0., fd_rate=0., Tsys=535.,
#                       CN0=45., ϕ=0., nADC=4, B=2.046e7,
#                       include_carrier=true, include_adc=true,
#                       include_noise=true, code_start_idx=1,
#                       include_databits=true)
#     sample_num = Int(f_s * t_length)
#     # Generate time vector
#     t = Array{Float64}(undef, sample_num)
#     @threads for i in 1:sample_num
#         @inbounds t[i] = (i-1)/f_s  # s
#     end
#     # Calculate number of data bits (use `ceil` to ensure that databits do not repeat)
#     db_code_length = Int64(ceil(L5_db_chipping_rate*t_length))
#     # Calculate code chipping rates with Doppler applied
#     # L5I
#     f_l5i_d = L5_chipping_rate*(1. + f_d/L5_freq)
#     f_l5i_dd = L5_chipping_rate*fd_rate/L5_freq
#     # Neuman sequence
#     f_nh_d = nh_chipping_rate*(1. + f_d/L5_freq)
#     f_nh_dd = nh_chipping_rate*fd_rate/L5_freq
#     # Calculate the L5Q and nh code phase offsets
#     l5q_init_code_phase = calcinitcodephase(L5_code_length,
#                                             f_l5i_d, f_l5i_dd,
#                                             f_s, code_start_idx)
#     nh_init_code_phase = calcinitcodephase(nh_code_length,
#                                            f_nh_d, f_nh_dd,
#                                            f_s, code_start_idx)
#     # Allocate space for signal
#     data = Array{Complex{Float64}}(undef, sample_num)
#     isreplica = false
#     noexp = false
#     # Generate random databit vector
#     databits = rand(0:1, db_code_length)
#     return L5QSignal(type, prn, f_s, t_length, f_if, f_d, fd_rate,
#                      Tsys, CN0, ϕ, nADC, B, code_start_idx,
#                      l5q_init_code_phase, nh_init_code_phase, t,
#                      data, include_carrier, include_adc,
#                      include_noise,
#                      f_l5q_d, f_l5q_dd,
#                      f_nh_d, f_nh_dd, sample_num,
#                      isreplica, noexp, L5_chipping_rate,
#                      L5_freq, L5_code_length, db_code_length, databits)
# end
#
#
# """
#     definesignal!(signal::L5ISignal;
#                   prn=signal.prn, f_d=signal.f_d,
#                   f_if=signal.f_if, fd_rate=signal.fd_rate,
#                   Tsys=signal.Tsys, CN0=signal.CN0,
#                   ϕ=signal.ϕ, nADC=signal.nADC,
#                   B=signal.B,
#                   include_carrier=signal.include_carrier,
#                   include_adc=signal.include_adc,
#                   include_noise=signal.include_noise,
#                   code_start_idx=signal.code_start_idx,
#                   isreplica=signal.isreplica,
#                   include_databits=signal.include_databits)

# Redefine properties of locally generated L5Q signal
# based off its type, PRN, etc. 

# Sampling rate (f_s) and signal length (t_length)
# are the only parameters that cannot be redefined.

# """
# function definesignal!(signal::L5ISignal;
#                        prn=signal.prn, f_d=signal.f_d,
#                        f_if=signal.f_if, fd_rate=signal.fd_rate,
#                        Tsys=signal.Tsys, CN0=signal.CN0,
#                        ϕ=signal.ϕ, nADC=signal.nADC,
#                        B=signal.B,
#                        include_carrier=signal.include_carrier,
#                        include_adc=signal.include_adc,
#                        include_noise=signal.include_noise,
#                        code_start_idx=signal.code_start_idx,
#                        isreplica=signal.isreplica,
#                        noexp=signal.noexp,
#                        include_databits=signal.include_databits)
#     ## Calculate code chipping rates with Doppler applied
#     # L5Q
#     f_l5i_d = L5_chipping_rate*(1. + f_d/L5_freq)
#     f_l5i_dd = L5_chipping_rate*fd_rate/L5_freq
#     # Neuman sequence
#     f_nh_d = nh_chipping_rate*(1. + f_d/L5_freq)
#     f_nh_dd = nh_chipping_rate*fd_rate/L5_freq
#     # Calculate the L5Q and nh code phase offsets
#     l5i_init_code_phase = calcinitcodephase(L5_code_length,
#                                             f_l5i_d, f_l5i_dd,
#                                             signal.f_s,
#                                             code_start_idx)
#     nh_init_code_phase = calcinitcodephase(nh_code_length,
#                                            f_nh_d, f_nh_dd,
#                                            signal.f_s,
#                                            code_start_idx)
#     # Store udated variables to "L5QSignal" struct
#     signal.prn = prn
#     signal.f_d = f_d
#     signal.f_if = f_if
#     signal.fd_rate = fd_rate
#     signal.Tsys = Tsys
#     signal.CN0 = CN0
#     signal.ϕ = ϕ
#     signal.nADC = nADC
#     signal.B = B
#     signal.include_carrier = include_carrier
#     signal.include_adc = include_adc
#     signal.include_noise = include_noise
#     signal.f_l5q_d = f_l5q_d
#     signal.f_l5q_dd = f_l5q_dd
#     signal.f_nh_d = f_nh_d
#     signal.f_nh_dd = f_nh_dd
#     signal.l5q_init_code_phase = l5q_init_code_phase
#     signal.nh_init_code_phase = nh_init_code_phase
#     signal.code_start_idx = code_start_idx
#     signal.isreplica = isreplica
#     signal.noexp = noexp
#     signal.include_databits = include_databits
#     return signal
# end
