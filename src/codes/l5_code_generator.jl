"""
    L5_code_length


The number of bits in the L5 signal for a 1ms segment = `10230`
"""
const L5_code_length = 10230  # chips


"""
    nh20_code_length


The number of bits in the Neuman sequence in a 1s segment = `20`
"""
const nh20_code_length = 20  # chips


"""
    nh10_code_length


The number of bits in the Neuman sequence in a 1s segment = `10`
"""
const nh10_code_length = 10  # chips


"""
    L5_chipping_rate


The chipping rate of the L5 signal = `10.23e6Hz`
"""
const L5_chipping_rate = 10.23e6  # Hz


"""
    nh_chipping_rate


The chipping rate of the Neuman sequence = `1000Hz`
"""
const nh_chipping_rate = 1000.  # Hz


"""
   L5_db_chipping_rate


The chipping rate of the databits for L5I = `100Hz`
"""
const L5_db_chipping_rate = 100.  # Hz


"""
    XA_XB_code_length


The length of the XA and XB codes used to construct the L5 signal = `8190`
"""
const XA_XB_code_length = 8190  # chips


"""
    XA_taps


The taps for the XA shift register.

Taps used = `[9, 10, 12, 13]`
"""
const XA_taps = [9, 10, 12, 13]


"""
    XB_taps


The taps for the XB shift register.

Taps used = `[1, 3, 4, 6, 7, 8, 12, 13]`
"""
const XB_taps = [1, 3, 4, 6, 7, 8, 12, 13]


"""
    XB_Q_intial_conditions

The initial conditions of the XB shift register per PRN for Q5 only. Used with
the dataless signal.
"""
const XB_Q_intial_conditions = Dict( 1 => [1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0],
                                     2 => [0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0],
                                     3 => [1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1],
                                     4 => [0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0],
                                     5 => [0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0],
                                     6 => [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1],
                                     7 => [1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1],
                                     8 => [0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0],
                                     9 => [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1],
                                    10 => [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0],
                                    11 => [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1],
                                    12 => [0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1],
                                    13 => [0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1],
                                    14 => [1, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                                    15 => [1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1],
                                    16 => [1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1],
                                    17 => [1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0],
                                    18 => [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0],
                                    19 => [0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1],
                                    20 => [1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1],
                                    21 => [0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
                                    22 => [0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0],
                                    23 => [1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1],
                                    24 => [0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1],
                                    25 => [0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1],
                                    26 => [0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0],
                                    27 => [1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0],
                                    28 => [1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0],
                                    29 => [0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0],
                                    30 => [1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1],
                                    31 => [0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1],
                                    32 => [1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0])


"""
    XB_I_intial_conditions


The initial conditions of the XB shift register per PRN for I5 only. Used with
the data signal.
"""
const XB_I_intial_conditions = Dict( 1 => [0, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0],
                                     2 => [1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1],
                                     3 => [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                                     4 => [1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0],
                                     5 => [1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1],
                                     6 => [0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0],
                                     7 => [1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1],
                                     8 => [1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0],
                                     9 => [1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1],
                                    10 => [0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0],
                                    11 => [0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0],
                                    12 => [1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1],
                                    13 => [0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0],
                                    14 => [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1],
                                    15 => [0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0],
                                    16 => [0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1],
                                    17 => [0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1],
                                    18 => [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0],
                                    19 => [1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1],
                                    20 => [0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1],
                                    21 => [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                                    22 => [1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1],
                                    23 => [1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0],
                                    24 => [1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0],
                                    25 => [1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1],
                                    26 => [1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0],
                                    27 => [0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0],
                                    28 => [0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0],
                                    29 => [0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1],
                                    30 => [1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],
                                    31 => [0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0],
                                    32 => [0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1])


"""
    nh20


A 20-bit Neuman-Hofman code XORed to the L5Q code. `nh20` is 20ms in length,
where every bit is XORed to each 1ms L5Q code sequence.
"""
const nh20 = [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0]


"""
    nh10


A 10-bit Neuman-Hofman code XORed to the L5I code. `nh10` is 10ms in length,
where every bit is XORed to each 1ms L5I code sequence.
"""
const nh10 = [0, 0, 0, 0, 1, 1, 0, 1, 0, 1]


"""
    xaregister()


Generates the XA register 8091 bit sequence. Restarts after the first 8090 bits,
where it repeats a second time until the 2040th state where it is reset again
for the start of the next 1ms.


Arguments: `None`


Returns:

- `xa_code::Vector`: vector containg the XA register result for 1ms period of
                     time
"""
function xaregister()
    xa_register = ones(Int, 13)
    xa_code = zeros(Int, L5_code_length)
    @inbounds for i in 1:XA_XB_code_length
        xa_code[i] = xa_register[end]
        xa_register[end] = xor(xa_register[XA_taps]...)
        xa_register = circshift(xa_register, 1)
    end
    xa_code[XA_XB_code_length+1:end] = xa_code[1:2040]
    return xa_code
end

"""
    xbregister(prn::Int, intial_conditions)


Generates the XB code sequence for 1ms for the specified PRN. Register starts
with initial conditions defined by `intial_conditions[prn]`.


Required Arguments:

- `prn::Int`: the pseudorandom noise (PRN) number, which is also related to the
              satellite vehicle
- `initial_conditions::Dict`: dictionary containing the initial states of the
                              register for each PRN


Returns:

- `xb_code::Vector`: vector containing the XB code sequence for a given PRN
                     and initial conditions
"""
function xbregister(prn::Int, intial_conditions)
    xb_register = deepcopy(intial_conditions[prn])
    xb_code = zeros(eltype(intial_conditions[prn]),
    	            L5_code_length)
    @inbounds for i in 1:L5_code_length
        xb_code[i] = xb_register[end]
        xb_register[end] = xor(xb_register[XB_taps]...)
        xb_register = circshift(xb_register, 1)
    end
    return xb_code
end


"""
    xa_code


Constant output of the XA register for the L5Q signal.
"""
const xa_code = xaregister()


"""
    gen_xb_codes()


Generates the output of the XB registers for all PRNs in `XB_intial_conditions.`


Arguments: `None`


Returns:

- `xb_codes::Dict`: dictionary containing the output of the XB register for all
                    PRN numbers
"""
    function gen_xb_codes()
	xb_codes = typeof(XB_Q_intial_conditions)()
	for prn in keys(XB_Q_intial_conditions)
		xb_codes[prn] = xbregister(prn, XB_Q_intial_conditions)
	end
	return xb_codes
end


# """
#     xb_codes

# Dictionary containing the output of the XB register for
# all PRNs.
# """
# const xb_codes = gen_xb_codes(XB_intial_conditions)


"""
    gen_l5q_codes()


Generates the L5Q codes for each PRN listed in `XB_Q_intial_conditions.`


Arguments: `None`


Returns:

- `l5q_codes::Dict`: dictionary containing L5Q codes for all PRN numbers
"""
function gen_l5q_codes()
	l5q_codes = typeof(XB_Q_intial_conditions)()
	for prn in keys(XB_Q_intial_conditions)
		l5q_codes[prn] = xor.(xa_code, xbregister(prn, XB_Q_intial_conditions))
	end
	return l5q_codes
end


"""
    gen_l5i_codes()


Generates the L5I codes for each PRN listed in `XB_I_intial_conditions.`


Arguments: `None`


Returns:

- `l5i_codes::Dict`: dictionary containing L5I codes for all PRN numbers
"""
function gen_l5i_codes()
    l5i_codes = typeof(XB_I_intial_conditions)()
    for prn in keys(XB_I_intial_conditions)
        l5i_codes[prn] = xor.(xa_code, xbregister(prn, XB_I_intial_conditions))
    end
    return l5i_codes
end


"""
    l5q_codes


Dictionary containing the L5Q codes for all PRNs.
"""
const l5q_codes = gen_l5q_codes()


"""
    l5i_codes


Dictionary containing the L5I codes for all PRNs.
"""
const l5i_codes = gen_l5i_codes()


"""
	define_l5_code_type(t_length=missing; channel="both",
						prns=collect(keys(l5i_codes)), sig_freq=L5_freq)


Generates a `SignalType` struct which defines the L5 codes that can be used
to define and generate L5 signals. The databits generated are random.


Optional Arguments:

- `t_length:Float64`: the length of the planned signal in seconds
    * this argument is used only to generate a `Vector` of databits so they do
      not repeat throughout the duration of the signal, once it is generated
    * `(default = missing)`
	* if given, databits are included in signal type
- `channel`: channel to include (either `"I"`, `"Q"`, or `"both"`)
	* `(default = "both")`
- `prns`: codes to define for signal type `(default = keys(l5i_codes))`
	* default is to define all PRN codes and assume that databits for each
	  are different
- `sig_freq`: carrier frequency in Hz `(default = L5_freq)`
- `databits`: optional databits `(default = missing)`, where random bits are used
- `B`: receiver bandwidth `(default = missing)`


Returns:

- `siganl_type::SignalType`: can be used in `definesignal` function to create
                             a signal, which can be used with `generatesignal!`
                             to generate the signal
"""
function define_l5_code_type(t_length=missing; channel="both",
	                         prns=collect(keys(l5i_codes)), sig_freq=L5_freq,
							 databits=missing, B=missing)
	if (channel == "I") || (channel == "both")
        if ~ismissing(t_length)
            # Generate random databits
			if ismissing(databits)
				databits = random_databits(L5_db_chipping_rate, t_length;
				                           prns=prns)
	        end
            # Include databits in I channel codes
            I_codes = definecodetype([copy_dictionary(l5i_codes, prns), nh10],
			                         [L5_chipping_rate, nh_chipping_rate];
                                     databits=[databits, L5_db_chipping_rate])
	    else
	        # Define I channel codes without databits
	        I_codes = definecodetype([copy_dictionary(l5i_codes, prns), nh10],
			                         [L5_chipping_rate, nh_chipping_rate])
	    end
		if channel == "I"
			signal_type = definesignaltype(I_codes, sig_freq, "I"; name="L5",
                                           B=B)
		else  # means that `channel` == "both"
			Q_codes = definecodetype([copy_dictionary(l5q_codes, prns), nh20],
									 [L5_chipping_rate, nh_chipping_rate])
			signal_type = definesignaltype(I_codes, Q_codes, sig_freq; name="L5",
                                           B=B)
		end
	elseif channel == "Q"
		Q_codes = definecodetype([copy_dictionary(l5q_codes, prns), nh20],
								 [L5_chipping_rate, nh_chipping_rate])
		signal_type = definesignaltype(Q_codes, sig_freq, "Q"; name="L5", B=B)
	else
		error("Invalid channel specifed.")
	end
    return signal_type
end
