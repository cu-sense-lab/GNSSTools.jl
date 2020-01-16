## L5Q constants
"""
L5 Frequency (Hz)
"""
const L5_freq = 1176.45e6  # Hz


"""
L5_code_length

The number of bits in the L5 signal
for a 1ms segment.
"""
const L5_code_length = 10230  # chips


"""
nh_code_length

The number of bits in the Neuman
sequence in a 1s segment.
"""
const nh_code_length = 20  # chips

"""
L5_chipping_rate

The chipping rate of the L5 signal
in Hz.
"""
const L5_chipping_rate = 10.23e6  # Hz


"""
nh_chipping_rate

The chipping rate of the Neuman
sequence in Hz
"""
const nh_chipping_rate = 1000.  # Hz


"""
XA_XB_code_length

The length of the XA and XB codes used to
construct the L5 signal.
"""
const XA_XB_code_length = 8190  # chips


"""
XA_taps

The taps for the XA shift register.
"""
const XA_taps = [9, 10, 12, 13]


"""
XB_taps

The taps for the XB shift register.
"""
const XB_taps = [1, 3, 4, 6, 7, 8, 12, 13]


"""
XB_intial_conditions

The initial conditions of the XB shift register
per PRN for Q5 only.
"""
const XB_intial_conditions = Dict( 1 => [1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0],
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
nh20

A 20-bit Neuman-Hofman code XORed to the L5 Q code.
nh20 is 20ms in length, where every bit is XORed to
each 1ms L5 Q code sequence.
"""
const nh20 = [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0]


"""
xaregister()

Generates the XA register 8091 bit sequence.
Restarts after the first 8090 bits, where it
repeats a second time until the 2040th state
where it is reset again for the start of the
next 1ms.
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
xbregister(prn::Integer)

Generates the XB code sequence for 1ms for the
specified PRN.
"""
function xbregister(prn::Integer)
    xb_register = deepcopy(XB_intial_conditions[prn])
    xb_code = zeros(eltype(XB_intial_conditions[prn]),
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

The output of the XA register for the L5Q signal.
"""
const xa_code = xaregister()


"""
gen_xb_codes()

Generates the output of the XB registers for all PRNs
in "XB_intial_conditions."
"""
function gen_xb_codes(XB_intial_conditions)
	xb_codes = typeof(XB_intial_conditions)()
	for prn in keys(XB_intial_conditions)
		xb_codes[prn] = xbregister(prn)
	end
	return xb_codes
end


# """
# xb_codes

# Dictionary containing the output of the XB register for
# all PRNs.
# """
# const xb_codes = gen_xb_codes(XB_intial_conditions)


"""
gen_l5q_codes()


Generates the L5Q codes for each PRN listed in
"XB_intial_conditions."
"""
function gen_l5q_codes(XB_intial_conditions)
	l5q_codes = typeof(XB_intial_conditions)()
	for prn in keys(XB_intial_conditions)
		l5q_codes[prn] = xor.(xa_code, xbregister(prn))
	end
	return l5q_codes
end


# """
# l5q_codes

# Dictionary containing the L5Q codes for all PRNs.
# """
# const l5q_codes = gen_l5q_codes(XB_intial_conditions)
