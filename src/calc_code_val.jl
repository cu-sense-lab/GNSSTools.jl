"""
    binary2OneOrNegativeOne(val) 


Converts an integer that is either 0 or 1 to -1 or 1. Values that are 0 are 
returned as -1 and values that are 1 are returned as 1.


Required Arguments:

- `val`: integer value that is either 0 or 1


Returns:

- a value that is either -1 or 1
"""
function binary2OneOrNegativeOne(val) 
    if val == 1
        return 1
    elseif val == 0
        return -1
    else
        error("Value must be either 0 or 1.")
    end
end


"""
    calc_code_val(signal::ReplicaSignal, t)


Calculates the value of the generic code with parameters defined by `signal` at
a time specified by `t` in seconds. Returns a Complex{Int} value that is either
-1 or 1 and a phase in `rad`. Use the `code_ϕ`, or carrier phase` result to
perform QPSK correctly.


Required Arguments:

- `signal::ReplicaSignal`: struct containing signal data and parameters
- `t::Float64`: current time in seconds


Returns:

- `(code_val, code_ϕ)::Tuple` where
    * `code_val::Complex{Int}`: unnormalized complex number with the values of
                                the I and Q channel codes being the real and
                                imaginary parts of this value, respectivly
    * `code_ϕ::Float64`: the carrier phase in `rad` corresponding to the values
                         of the I and Q channel codes


The carrier phase value for a given set of I and Q channel code values at a
given time `t` is given by:

```math
    code_\\phi = \\tan\\left(\\frac{Q}{I}\\right)^{-1}
```

where for codes on **both** the I and Q channel the phase is:

| I Code Value  | Q Code Value  | Phase (degrees)  |
|--------------:|--------------:|-----------------:|
| 1             | 1             | 45               |
| -1            | 1             | 135              |
| -1            | -1            | 225              |
| 1             | -1            | 315              |

and for codes on **only** the I channel:

| I Code Value  | Q Code Value  | Phase (degrees)  |
|--------------:|--------------:|-----------------:|
| 1             | 0             | 0                |
| -1            | 0             | 180              |

and for codes on **only** the Q channel:

| I Code Value  | Q Code Value  | Phase (degrees) |
|--------------:|--------------:|----------------:|
| 0             | 1             | 90              |
| 0             | -1            | 270             |
"""
function calc_code_val(signal::ReplicaSignal, t)
    prn = signal.prn
    # I channel
    # If I channel codes exist, Multiply the value of each layer of code at a given
    # `t` with the value of `I_val`
    if signal.signal_type.include_I
        # I_val = 0
        I_val = 1
        I_codes = signal.signal_type.I_codes
        # Determine whether databits are included in I channel and whether
        # `include_databits_I` flag is true. If so, databits are multiplied by
        # the value of `I_val`. If this flag is set to false or if
        # `signal.isreplica` is set to true, no databits are included.
        if I_codes.databits
            if ~signal.isreplica && signal.include_databits_I
                code_num = I_codes.code_num
            else
                code_num = I_codes.code_num - 1
            end
        else
            code_num = I_codes.code_num
        end
        # Loop through I channel codes and multiply current value at time `t` to
        # value of `I_val`
        for i in 1:code_num
            code_idx = calccodeidx(signal.init_code_phases_I[i],
                                   signal.f_code_d_I[i], signal.f_code_dd_I[i],
                                   t, I_codes.code_lengths[i])
            # I_val = xor(I_val, I_codes.codes[i][prn][code_idx])
            I_val *= binary2OneOrNegativeOne(I_codes.codes[i][prn][code_idx])
        end
        # Change from 0:1 state to -1:1 state. This prevensts divide over zero
        # error in `atan` function at the end of this function.
        # I_val = binary2OneOrNegativeOne(I_val)
    else
        I_val = 0
    end
    # Q channel
    # If Q channel codes exist, Multiply the value of each layer of code at a given
    # `t` with the value of `Q_val`
    if signal.signal_type.include_Q
        # Q_val = 0
        Q_val = 1
        Q_codes = signal.signal_type.Q_codes
        # Determine whether databits are included in Q channel and whether
        # `include_databits_Q` flag is true. If so, databits are multiplied by
        # the value of `Q_val`. If this flag is set to false or if
        # `signal.isreplica` is set to true, no databits are included.
        if Q_codes.databits
            if ~signal.isreplica && signal.include_databits_Q
                code_num = Q_codes.code_num
            else
                code_num = Q_codes.code_num - 1
            end
        else
            code_num = Q_codes.code_num
        end
        # Loop through Q channel codes and multiply current value at time `t` to
        # value of `Q_val`
        for i in 1:code_num
            code_idx = calccodeidx(signal.init_code_phases_Q[i],
                                   signal.f_code_d_Q[i], signal.f_code_dd_Q[i],
                                   t, Q_codes.code_lengths[i])
            # Q_val = xor(Q_val, Q_codes.codes[i][prn][code_idx])
            Q_val *= binary2OneOrNegativeOne(Q_codes.codes[i][prn][code_idx])
        end
        # Change from 0:1 state to -1:1 state. This prevensts divide over zero
        # error in `atan` function at the end of this function.
        # Q_val = binary2OneOrNegativeOne(Q_val)
    else
        Q_val = 0
    end
    # Calculate the phase of the carrier based off the values of `Q_val` and
    # `I_val`. This simulates a QPSK signal.
    code_ϕ = atan(Q_val, I_val)
    # Construct a `Complex{Int}` that has the value of the I channel codes in
    # `real(code_val)` and the value of the Q channel codes in `iamg(code_val)`.
    code_val = complex(I_val, Q_val)
    return (code_val, code_ϕ)
end


"""
    calc_code_val(signal::ReplicaSignal, t, code_chips_I, code_chips_Q)


Calculates the value of the generic code with parameters defined by `signal` at
a time specified by `t` in seconds. Returns a Complex{Int} value that is either
-1 or 1 and a phase in `rad`. Use the `code_ϕ`, or carrier phase` result to
perform QPSK correctly.

This method accepts arrays `code_chips_I` and `code_chips_Q`, which are arrays
containing interpolated functions that calculate the current code chip as a
function of `t` in seconds. This method is used when the user supplies
`generatesignal!` with `doppler_curve` and `doppler_t` arrays, where
`doppler_curve` is assumed to be a non-linear Doppler frequency curve.


Required Arguments:

- `signal::ReplicaSignal`: struct containing signal data and parameters
- `t::Float64`: current time in seconds
- `code_chips_I`: N element array containing interpolated functions
                  that compute the current code chip for a given `t`
                  for each code on the I channel
- `code_chips_Q`: M element array containing interpolated functions
                  that compute the current code chip for a given `t`
                  for each code on the Q channel


Returns:

- `(code_val, code_ϕ)::Tuple` where
    * `code_val::Complex{Int}`: unnormalized complex number with the values of
                                the I and Q channel codes being the real and
                                imaginary parts of this value, respectivly
    * `code_ϕ::Float64`: the carrier phase in `rad` corresponding to the values
                         of the I and Q channel codes


The carrier phase value for a given set of I and Q channel code values at a
given time `t` is given by:

```math
    code_\\phi = \\tan\\left(\\frac{Q}{I}\\right)^{-1}
```

where for codes on **both** the I and Q channel the phase is:

| I Code Value  | Q Code Value  | Phase (degrees)  |
|--------------:|--------------:|-----------------:|
| 1             | 1             | 45               |
| -1            | 1             | 135              |
| -1            | -1            | 225              |
| 1             | -1            | 315              |

and for codes on **only** the I channel:

| I Code Value  | Q Code Value  | Phase (degrees)  |
|--------------:|--------------:|-----------------:|
| 1             | 0             | 0                |
| -1            | 0             | 180              |

and for codes on **only** the Q channel:

| I Code Value  | Q Code Value  | Phase (degrees) |
|--------------:|--------------:|----------------:|
| 0             | 1             | 90              |
| 0             | -1            | 270             |
"""
function calc_code_val(signal::ReplicaSignal, t, code_chips_I, code_chips_Q)
    prn = signal.prn
    # I channel
    # If I channel codes exist, Multiply the value of each layer of code at a
    # given `t` with the value of `I_val`
    if signal.signal_type.include_I
        # I_val = 0
        I_val = 1
        I_codes = signal.signal_type.I_codes
        # Determine whether databits are included in I channel and whether
        # `include_databits_I` flag is true. If so, databits are multiplied by
        # the value of `I_val`. If this flag is set to false or if
        # `signal.isreplica` is set to true, no databits are included.
        if I_codes.databits
            if ~signal.isreplica && signal.include_databits_I
                code_num = I_codes.code_num
            else
                code_num = I_codes.code_num - 1
            end
        else
            code_num = I_codes.code_num
        end
        # Loop through I channel codes and multiply current value at time `t` to
        # value of `I_val`
        for i in 1:code_num
            code_chip = code_chips_I[i](t)
            code_idx = Int(floor(code_chip%signal.signal_type.I_codes.code_lengths[i])) + 1
            # I_val = xor(I_val, I_codes.codes[i][prn][code_idx])
            I_val *= binary2OneOrNegativeOne(I_codes.codes[i][prn][code_idx])
        end
        # Change from 0:1 state to -1:1 state. This prevensts divide over zero
        # error in `atan` function at the end of this function.
        # I_val = binary2OneOrNegativeOne(I_val)
    else
        I_val = 0
    end
    # Q channel
    # If Q channel codes exist, Multiply the value of each layer of code at a 
    # given `t` with the value of `Q_val`
    if signal.signal_type.include_Q
        # Q_val = 0
        Q_val = 1
        Q_codes = signal.signal_type.Q_codes
        # Determine whether databits are included in Q channel and whether
        # `include_databits_Q` flag is true. If so, databits are multiplied by
        # the value of `Q_val`. If this flag is set to false or if
        # `signal.isreplica` is set to true, no databits are included.
        if Q_codes.databits
            if ~signal.isreplica && signal.include_databits_Q
                code_num = Q_codes.code_num
            else
                code_num = Q_codes.code_num - 1
            end
        else
            code_num = Q_codes.code_num
        end
        # Loop through Q channel codes and multiply current value at time `t` to
        # value of `Q_val`
        for i in 1:code_num
            code_chip = code_chips_Q[i](t)
            code_idx = Int(floor(code_chip%signal.signal_type.Q_codes.code_lengths[i])) + 1
            # Q_val = xor(Q_val, Q_codes.codes[i][prn][code_idx])
            Q_val *= binary2OneOrNegativeOne(Q_codes.codes[i][prn][code_idx])
        end
        # Change from 0:1 state to -1:1 state. This prevensts divide over zero
        # Q_val = binary2OneOrNegativeOne(Q_val)
    else
        I_val = 0
    end
    # Calculate the phase of the carrier based off the values of `Q_val` and
    # `I_val`. This simulates a QPSK signal.
    code_ϕ = atan(Q_val, I_val)
    # Construct a `Complex{Int}` that has the value of the I channel codes in
    # `real(code_val)` and the value of the Q channel codes in `iamg(code_val)`.
    code_val = complex(I_val, Q_val)
    return (code_val, code_ϕ)
end
