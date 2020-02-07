"""
    Klobuchar(α₀, α₁, α₁, α₃, β₀, β₁, β₂, β₃)

Struct holding the ionospheric correction terms
for the Klobuchar Ionospheric Model.
"""
struct Klobuchar
    α₀::Float64
    α₁::Float64
    α₂::Float64
    α₃::Float64
    β₀::Float64
    β₁::Float64
    β₂::Float64
    β₃::Float64
end


"""
    BRDC()

Struct holding GPS Broadcast Ephemeris parameters and
ionospheric correction terms.
"""
struct BRDC{T}
    prn::Int64
    iono_parms::T
    Af₀::Float64
    Af₁::Float64
    Af₂::Float64
    IODE::Float64
    Crs::Float64
    Crc::Float64
    Cus::Float64
    Cuc::Float64
    Cis::Float64
    Cic::Float64
    ṅ::Float64
    M₀::Float64
    e::Float64
    a::Float64
    Toe::Float64
    Toc::Float64
    Ω::Float64
    i::Float64
    di::Float64
    dα::Float64
    ω::Float64
    gps_week::Float64
    health::Float64
end


"""
    replaceD(str)

Replaces a `D` with an `E` from values in RINEX files
so they can be parsed using `parse`.
"""
function replaceD(str)
    str = replace(str, "D" => "E")
    return str
end


"""
    loadrinex(file)

Loads a RINEX broadcast ephemeris file and stores into
dictionary where the keys are the PRN numbers.

Julia implementation based off the MATLAB implementation,
writen by Ben K. Bradley (2009-07-19).

Returns dictionary containing `BRDC` structs for each PRN.
The dictionary keys are the PRN values.
"""
function loadrinex(file)
    f = open(file, "r")
    # Read header
    endofheader = false
    alphas = missing
    betas = missing
    fline = missing
    while ~endofheader
        fline = readline(f)
        # Check if this is the last line of the header
        endofheader = occursin("END OF HEADER", fline)
        # Check if this is the line with the ionospheric correction terms, αᵢ
        if occursin("ION ALPHA", fline)
            alphas = split(fline)
            alphas = [parse(Float64, replaceD(alphas[1])),
                      parse(Float64, replaceD(alphas[2])),
                      parse(Float64, replaceD(alphas[3])),
                      parse(Float64, replaceD(alphas[4]))]
        end
        # Check if this is the line with the ionospheric correction terms, βᵢ
        if occursin("ION BETA", fline)
            betas = split(fline)
            betas = [parse(Float64, replaceD(betas[1])),
                     parse(Float64, replaceD(betas[2])),
                     parse(Float64, replaceD(betas[3])),
                     parse(Float64, replaceD(betas[4]))]
        end
    end
    iono_parms = Klobuchar(alphas[1], alphas[2], alphas[3], alphas[4],
                           betas[1], betas[2], betas[3], betas[4])
    BRDCs = Dict()
    fline = readline(f)  # Read line
    while fline != ""
        """
            First Line of PRN Navigation Data
        """
        prn = parse(Int64, fline[1:2], base=10)       # PRN
        # Read clock information
        Af₀ = parse(Float64, replaceD(fline[23:41]))  # clock bias (s)
        Af₁ = parse(Float64, replaceD(fline[42:60]))  # clock drift (s/s)
        Af₂ = parse(Float64, replaceD(fline[61:79]))  # clock drift rate (s/s²)
        """
            Second Line of PRN Navigation Data
        """
        fline = readline(f)                           # Read line
        IODE = parse(Float64, replaceD(fline[4:22]))  # Issue of data (ephemeris)
        Crs = parse(Float64, replaceD(fline[23:41]))  # Amplitude of the Sine Harmonic Correction
                                                      # Apply to orbital radius
        ṅ = parse(Float64, replaceD(fline[42:60]))    # Mean motion difference from computed value (rad/s)
        M₀ = parse(Float64, replaceD(fline[42:60]))   # Mean anomaly at reference time (rad)
        """
            Third Line of PRN Navigation Data
        """
        fline = readline(f)                          # Read line
        Cuc = parse(Float64, replaceD(fline[4:22]))  # Amplitude of the Cosine Harmonic Correction
                                                     # Apply to Argument of Latitude
        e = parse(Float64, replaceD(fline[23:41]))   # Eccentricity
        Cus = parse(Float64, replaceD(fline[42:60])) # Amplitude of the Sine Harmonic Correction
                                                     # Apply to Argument of Latitude
        a = parse(Float64, replaceD(fline[61:79]))^2 # Semi-major axis (meters)
        """
            Fourth Line of PRN Navigation Data
        """
        fline = readline(f)                          # Read line
        Toe = parse(Float64, replaceD(fline[4:22]))  # Referemce time of ephemeris (seconds into GPS week)
        Cic = parse(Float64, replaceD(fline[23:41])) # Amplitude of Cosine Harmonic Correction
                                                     # Apply to th Angle of Inclination (i)
        Ω = parse(Float64, replaceD(fline[42:60]))   # Longitude of Ascending Node of orbit plane
                                                     # at weekly epoch (rad)
        Cis = parse(Float64, replaceD(fline[61:79])) # Amplitude of Sine Harmonic Correction
                                                     # Apply to th Angle of Inclination (i)
        """
            Fifth Line of PRN Navigation Data
        """
        fline = readline(f)                          # Read line
        i = parse(Float64, replaceD(fline[4:22]))    # Inclination angle at reference time (rad)
        Crc = parse(Float64, replaceD(fline[23:41])) # Amplitude of the Cosine Harmonic Correction
                                                     # Use with orbt radius
        ω = parse(Float64, replaceD(fline[42:60]))   # Argument of perigee (rad)
        dα = parse(Float64, replaceD(fline[61:79]))  # Rate of change of Right Ascension (rad/s)
        """
            Sixth Line of PRN Navigation Data
        """
        fline = readline(f)                          # Read line
        di = parse(Float64, replaceD(fline[4:22]))   # Rate of change inclination angle (rad/s)
        week = parse(Float64, replaceD(fline[42:60]))# GPS week number (use with Toe)
        if week < 1024
            week += 1024
        end
        """
            Seventh Line of PRN Navigation Data
        """
        fline = readline(f)                          # Read line
        health = parse(Float64, replaceD(fline[23:41])) # Satellite health (0.00 = usable)
        """
            Eigtht Line of PRN Navigation Data
        """
        fline = readline(f)                          # Read line
        Toc = Toe                                    # Time of clock
        # Save to `BRDC` struct and place into dictionary, `BRDCs`
        BRDCs[prn] = BRDC(prn, iono_parms, Af₀, Af₁, Af₂, IODE, Crs, Crc,
                          Cus, Cuc, Cis, Cic, ṅ, M₀, e, a, Toe, Toc, Ω, i,
                          di, dα, ω, week, health)
        fline = readline(f)                           # Read line
    end
    return BRDCs
end