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

Returns the IF and sampling frequency for L1 frequency range file.
"""
function g1b1()
    f_s = 25e6  # samples-per-second
    f_center = 1.57e9  # Hz
    f_if = abs(L1_freq-f_center)
    sigtype = Val(:l1ca)
    return (f_s, f_if, f_center, L1_freq, sigtype)
end


"""
    g5()

Returns the IF and sampling frequency for L5 frequency range file.
"""
function g5()
    f_s = 25e6  # samples-per-second
    f_center = 1.17645e9  # Hz
    f_if = abs(L5_freq-f_center)
    sigtype = Val(:l5q)
    return (f_s, f_if, f_center, L5_freq, sigtype)
end



"""
    parseGPSData(data_file)

Parse the GPS data file from AGI.
Dictionary stored as `Dict(PRN, Dict(SV, Dict(other_fields...)))`.
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

Download GPS SV data from AGI and parse.
Returns a dictionary with historical SVs
for each PRN.
"""
function getCurrentGPSNORADIDS(GPS_data_file_url)
    directory = string(homedir(), "/.GNSSTools")
    if ~isdir(directory)
        mkpath(directory)
    end
    ids_file = string(directory, "/GPSData.txt")
    if ~isfile(ids_file)
        run(`curl -o $(ids_file) $(GPS_data_file_url)`);
    end
    GPSData = parseGPSData(ids_file)
    return GPSData
end


const GPS_data_file_url = "ftp://ftp.agi.com/pub/Catalog/Almanacs/SEM/GPSData.txt"


const GPSData = getCurrentGPSNORADIDS(GPS_data_file_url)
