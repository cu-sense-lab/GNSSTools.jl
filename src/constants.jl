"""
    Speed of light (m/s)

c = 299,792,458 m/s
"""
const c = 299792458  # m/s


"""
    Boltzman constant

k = 1.38×10⁻²³ J/K
"""
const k = 1.380649e-23  # J/K


"""
    Rₑ (Volumetric Mean Radius of the Earth)

Rₑ = 6371.000e3  # meters

From [NASA Earth Fact Sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
"""
const Rₑ = 6371.000e3  # meters


"""
    L1 frequency (Hz)
"""
const L1_freq = 1575.42e6  # Hz


"""
    L5 Frequency (Hz)
"""
const L5_freq = 1176.45e6  # Hz


"""
    g1b1()

Returns the IF and sampling frequency for L1 frequency range.
"""
function g1b1()
    f_s = 25e6  # samples-per-second
    f_center = 1.57e9  # Hz
    f_if = abs(L1_freq-f_center)
    return (f_s, f_if, f_center)
end


"""
    g5()

Returns the IF and sampling frequency for L5 frequency range.
"""
function g5()
    f_s = 25e6  # samples-per-second
    f_center = 1.17645e9  # Hz
    f_if = abs(L5_freq-f_center)
    return (f_s, f_if, f_center)
end
