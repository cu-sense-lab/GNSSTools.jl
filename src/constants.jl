"""
    c


Speed of light = `299,792,458m/s`
"""
const c = 299792458  # m/s


"""
    k


Boltzman constant = `1.38×10⁻²³ J/K`
"""
const k = 1.380649e-23  # J/K


"""
    Rₑ


Rₑ (Volumetric Mean Radius of the Earth) = `6371.000e3m`

From [NASA Earth Fact Sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html)
"""
const Rₑ = 6371.000e3  # meters


"""
    L1_freq


L1 carrier frequency = `1575.42e6Hz`
"""
const L1_freq = 1575.42e6  # Hz


"""
    L2_freq


L2 carrier frequency = `1227.60e6Hz`
"""
const L2_freq = 1227.60e6  # Hz


"""
    L5_freq


L5 Frequency `1176.45e6Hz`
"""
const L5_freq = 1176.45e6  # Hz


"""
    g1b1()


Data whose filenames contain the `g1b1` string are collected at 25 Msps,
with the center frequency set to 1.57GHz which may contain the following
signals:

- GPS L1
- Galileo E1bc
- Beidou B1i


Returns:

- `(f_s, f_if, f_center, L1_freq, sigtype)::Tuple`: where
    * `f_s::Float64`: sampling rate of data `(25MHz)`
    * `f_if::Float64`: absolute value of the difference between the signal
                       carrier frequency and receiver center frequency `(5.42MHz)`
    * `f_center::Float64`: receiver center frequency `(1.57GHz)`
    * `L1_freq::Float64`: GPS L1 carrier frequency `(1.57542GHz)`
    * `sigtype::Val{:l1ca}`: `[DEPRICATED]` used in previous implementation of
                             GNSSTools and will be removed soon
"""
function g1b1()
    f_s = 25e6  # Hz
    f_center = 1.57e9  # Hz
    f_if = abs(L1_freq-f_center)  # Hz
    sigtype = Val(:l1ca)
    return (f_s, f_if, f_center, L1_freq, sigtype)
end


"""
    g5()


Data whose filenames contain the `g5` string are collected at 25 Msps,
with the center frequency set to 1.17645GHz which may contain the following
signals:

- GPS L5
- Galileo E5a


Returns:

- `(f_s, f_if, f_center, L5_freq, sigtype)::Tuple`: where
    * `f_s::Float64`: sampling rate of data `(25MHz)`
    * `f_if::Float64`: absolute value of the difference between the signal
                       carrier frequency and receiver center frequency `(0Hz)`
    * `f_center::Float64`: receiver center frequency `(1.17645GHz)`
    * `L1_freq::Float64`: GPS L1 carrier frequency `(1.17645GHz)`
    * `sigtype::Val{:l5q}`: `[DEPRICATED]` used in previous implementation of
                            GNSSTools and will be removed soon
"""
function g5()
    f_s = 25e6  # Hz
    f_center = 1.17645e9  # Hz
    f_if = abs(L5_freq-f_center)  # Hz
    sigtype = Val(:l5q)
    return (f_s, f_if, f_center, L5_freq, sigtype)
end


"""
    WGS84


Struct type that hols the WGS84 constants.


Fields:

- `a::Float64`: semi-major axis of Earth `(6378137.0m)`
- `f::Float64`: flattening factor `(1/298.257223563)`
- `ω::Float64`: Earth's angular velocity `(7.2921151467e-5rad/s)`
- `μ::Float64`: also known as `GM`, is Earth's gravitational constant
                `(3.986005e14m³/s²)`
"""
struct WGS84
    a::Float64
    f::Float64
    ω::Float64
    μ::Float64
end


"""
    initWGS84()


Initializes the WGS84 constants and stores them in a `WGS84` struct.


Returns:

- `WGS84` struct
"""
function initWGS84()
    return WGS84(6378137.0,  # meters
                 1/298.257223563,
                 7.2921151467e-5,  # rad/s
                 3.986005e14)  # m³/s²
end


"""
    wgs84


Initialized constant struct that contains the WGS84 constants.
"""
const wgs84 = initWGS84()


"""
    parseGPSData(data_file)


Parse the GPS data file from AGI. See function `getCurrentGPSNORADIDS` below.
"""
function parseGPSData(data_file)
    file = open(data_file, "r")
    lines = readlines(file)
    close(file)
    N = size(lines)[1]
    GPSData = Dict{Int,Dict{Int,Dict{String,Any}}}()
    for i in 23:N
        line = lines[i]
        if length(line) > 1
            prn = parse(Int, line[3:4])
            sv = parse(Int, line[7:8])
            satnum = parse(Int, line[11:15])
            block = line[18:20]
            start_week = parse(Int, line[27:30])
            start_toa = parse(Int, line[32:37])
            start_year = parse(Int, line[39:42])
            start_month = parse(Int, line[44:45])
            start_day = parse(Int, line[47:48])
            start_date = (start_year, start_month, start_day)
            start_julian_date = DatetoJD(start_date..., 0, 0, 0)
            active = occursin("Active", line)
            notes = line[77:end]
            if active
                end_week = missing
                end_toa = missing
                end_date = missing
                end_julian_date = missing
            else
                end_week = parse(Int, line[53:56])
                end_toa = parse(Int, line[58:63])
                end_year = parse(Int, line[65:68])
                end_month = parse(Int, line[70:71])
                end_day = parse(Int, line[73:74])
                end_date = (end_year, end_month, end_day)
                end_julian_date = DatetoJD(end_date..., 0, 0, 0)
            end
            sv_info = Dict("satnum"=>satnum,
                           "block"=>block,
                           "start_week"=>start_week,
                           "start_toa"=>start_toa,
                           "start_date"=>start_date,
                           "start_julian_date"=>start_julian_date,
                           "end_week"=>end_week,
                           "end_toa"=>end_toa,
                           "end_date"=>end_date,
                           "end_julian_date"=>end_julian_date,
                           "active"=>active,
                           "status"=>notes)
            try
                GPSData[prn][sv] = sv_info
            catch
                GPSData[prn] = Dict(sv=>sv_info)
            end
        end
    end
    return GPSData
end


"""
    getCurrentGPSNORADIDS()

Download GPS SV data from AGI and parse. Returns a dictionary with historical
SVs for each PRN.


Required Arguments:

- `file_url::String`:
"""
function getCurrentGPSNORADIDS(file_url)
    directory = string(homedir(), "/.GNSSTools")
    if ~isdir(directory)
        mkpath(directory)
    end
    ids_file = string(directory, "/GPSData.txt")
    if ~isfile(ids_file)
        run(`curl -o $(ids_file) $(file_url)`);
    end
    GPSData = parseGPSData(ids_file)
    return GPSData
end


const GPS_data_file_url = "ftp://ftp.agi.com/pub/Catalog/Almanacs/SEM/GPSData.txt"


const GPSData = getCurrentGPSNORADIDS(GPS_data_file_url)
