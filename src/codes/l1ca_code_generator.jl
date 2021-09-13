"""
    l1ca_chipping_rate


L1 C/A Code chipping rate = `1.023e6Hz`
"""
const l1ca_chipping_rate = 1.023e6  # Hz


"""
    l1ca_chip_length


L1 C/A code chip length in 1ms = `1023`
"""
const l1ca_code_length = 1023  # chips


"""
    l1ca_db_chipping_rate


The chipping rate of the data bits for the L1 C/A code = `50Hz`
"""
const l1ca_db_chipping_rate = 50.  # Hz


"""
    l1ca_taps


Dictionary containing the phase selector taps (columns 1 and 2) and the code
chip delay amount in chips (column 3). The dictionary keys represent the PRNs.
"""
const l1ca_taps = Dict( 1 => (2,  6,   5),
                        2 => (3,  7,   6),
                        3 => (4,  8,   7),
                        4 => (5,  9,   8),
                        5 => (1,  9,  17),
                        6 => (2, 10,  18),
                        7 => (1,  8, 139),
                        8 => (2,  9, 140),
                        9 => (3, 10, 141),
                       10 => (2,  3, 251),
                       11 => (3,  4, 252),
                       12 => (5,  6, 254),
                       13 => (6,  7, 255),
                       14 => (7,  8, 256),
                       15 => (8,  9, 257),
                       16 => (9, 10, 258),
                       17 => (1,  4, 469),
                       18 => (2,  5, 470),
                       19 => (3,  6, 471),
                       20 => (4,  7, 472),
                       21 => (5,  8, 473),
                       22 => (6,  9, 474),
                       23 => (1,  3, 509),
                       24 => (4,  6, 512),
                       25 => (5,  7, 513),
                       26 => (6,  8, 514),
                       27 => (7,  9, 515),
                       28 => (8, 10, 516),
                       29 => (1,  6, 859),
                       30 => (2,  7, 860),
                       31 => (3,  8, 861),
                       32 => (4,  9, 862))


"""
    ca_code_chips(prn::Int)


Generates the CA Code chip-by-chip sequence for a given GPS satellite PRN
number.


Required Arguments:

- `prn::Int`: the pseudorandom noise (PRN) number, which is also related to the
              satellite vehicle


Returns:

- `ca_code::Vector{Int}`: a 1023-element array of C/A code chips for the PRN
                          number argument, `prn`
"""
function ca_code_chips(prn::Int)
    # Initialize G1 and G2 registers to 1's
    G1 = ones(Int64, 10)
    G2 = ones(Int64, 10)
    # Get phase selector taps for specific PRN
    phase_sel_taps = l1ca_taps[prn][1:2]
    ca_code = Array{Int64}(undef, l1ca_code_length)
    @inbounds for chip in 1:l1ca_code_length
        G1_output = G1[10]
        G1_epoch = xor(G1_output, G1[3])
        # Calculate phase selector output
        phase_sel_out = xor(G2[phase_sel_taps[1]],
                            G2[phase_sel_taps[2]])
        # XOR phase selector and G1 outputs to yield
        # C/A code value
        ca_code[chip] = xor(G1_output, phase_sel_out)
        # Calculate G2 epoch by XORing the G2 taps
        # 2, 3, 6, 8, 9, and 10
        G2_epoch = xor(G2[2], G2[3], G2[6], G2[8],
                       G2[9], G2[10])
        # Shift registers forward by one chip
        @inbounds for i in 10:-1:2
            G1[i] = G1[i-1]
            G2[i] = G2[i-1]
        end
        # Save G1 and G2 epoch values to beginning of
        # shifted G1 and G2 arrays
        G1[1] = G1_epoch
        G2[1] = G2_epoch
    end
    return ca_code
end


"""
    gen_l1ca_codes()


Generates a dictionary containing all the L1 C/A codes for PRNS 1 through 32.
The PRN number is the dictionary key.


Arguments: `None`


Returns:

- `l1ca_codes::Dict`: dictionary containing a 1023-element array of C/A code
                      chips for each PRN
"""
function gen_l1ca_codes()
    l1ca_codes = Dict{Int64, Array{Int64,1}}()
    for prn in keys(l1ca_taps)
        l1ca_codes[prn] = ca_code_chips(prn)
    end
    return l1ca_codes
end


"""
    l1ca_codes


Constant dictionary containing the L1 C/A codes. The keys are the PRN numbers.
"""
const l1ca_codes = gen_l1ca_codes()


"""
    define_l1ca_code_type(t_length=missing; prns=collect(keys(l1ca_codes)),
                          sig_freq=L1_freq)


Generates a `SignalType` struct which defines the L1 C/A code that can be used
to define and generate L1 C/A signals. The databits generated are random.


Optional Arguments:

- `t_length:Float64`: the length of the planned signal in seconds
    * this argument is used only to generate a `Vector` of databits so they do
      not repeat throughout the duration of the signal, once it is generated
    * `(default = missing)`
    * if given, databits are included in signal type
- `prns`: codes to define for signal type `(default = keys(l1ca_codes))`
	* default is to define all PRN codes and assume that databits for each
	  are different
- `sig_freq`: carrier frequency in Hz `(default = L1_freq)`
- `databits`: optional databits `(default = missing)`, where random bits are used
- `B`: receiver bandwidth `(default = missing)`


Returns:

- `siganl_type::SignalType`: can be used in `definesignal` function to create
                             a signal, which can be used with `generatesignal!`
                             to generate the signal
"""
function define_l1ca_code_type(t_length=missing; prns=collect(keys(l1ca_codes)),
                               sig_freq=L1_freq, databits=missing, B=missing)
    if ~ismissing(t_length)
        # Generate random databits
        if ismissing(databits)
            databits = random_databits(l1ca_db_chipping_rate, t_length; prns=prns)
        end
        # Include databits in I channel codes
        I_codes = definecodetype(copy_dictionary(l1ca_codes, prns),
                                 l1ca_chipping_rate;
                                 databits=[databits, l1ca_db_chipping_rate])
    else
        # Define I channel codes without databits
        I_codes = definecodetype(copy_dictionary(l1ca_codes, prns),
                                 l1ca_chipping_rate)
    end
    signal_type = definesignaltype(I_codes, sig_freq, "I"; name="L1 C/A", B=B)
    return signal_type
end
