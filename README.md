# GNSSTools

A general Julia package for working with GNSS data. Support course acquisition for L5Q signals. New will signals will be added on an as needed basis. Local signal generation and course acquisition methods support multithreading. See instructions for installing this package.

# Installation

## Installing GNSSTools

Begin by installing PyPlot, which the Python Matplotlib plotting library, ported to Julia. In the Julia terminal, type the following:

```Julia
ENV["PYTHON"]=""

```

Then, while still in Julia, type `]` to enter the package manager and copy/paste the following to setup PyPlot.


```Julia
add PyPlot
build PyPlot

```


After PyPlot is built, you can install GNSSTools using the following. Copy and paste this line while in the Julia package manager.


```julia
add https://github.com/cu-sense-lab/GNSSTools.jl.git
```

Julia will automatically install the neccessary dependencies and GNSSTools. Once complete, you can return the Julia REPL by hitting the `backspace` key. To update GNSSTools if newer versions are released, enter the package manager again by typing `]` and use `update GNSSTools`.


## TODOS

- [x] Simulate Starlink satellites (See [LabVideo Gitlab Repo](http://192.168.3.66/bilardis/LabVideo))
- [x] partial sample code offset
    * up-sample signal by N times
        + for each sample, N sub-samples are calculated and summed together
        + will not store N sub-samples; will instead calculate them in place and store only the sum of them
- [x] Use h-parameters to produce oscillator noise instead of the Voss-McCartney method
- [ ] Multi-PRN simulation
    * each signal from a given PRN is a modulated carrier of specific C/N₀
    * add multiple carriers together into a single raw IQ vector
- [x] Test LEO sat acquisition methods
- [ ] Generate L1C codes
- [x] Monte Carlo simulation for processing method evaluation
    * static C/N₀ and constellation
    * all other parameters are variable
        + thermal and phase noise are generated for each simulation run
        + initial code phase and carrier phase are picked from a uniform distribution at each iteration
- [ ] User dynamics due to user motion
    * allow user to provide a vector of user positions along with accompanying time vector
    * this vector will be used along with satellite information to calculate Doppler frequency curve
    * Doppler curve will be used for non-linear Doppler simulation
- [ ] Allow seed input for thermal and phase noise sources
- [ ] Effects due to ionosphere
    * not needed for X-band frequency, since ionospheric effects are negligible at that frequency
- [ ] Specify power ratio between I and Q channels
    * usually split 50% on I and Q
    * depending on code and band, it may be 75% on I and 25% on Q
- [ ] Track both I and Q channels simultaneously
    * may be able to do this by generating a replica with both I and Q codes
    * during tracking, track each channel separately
    * assign separate early, prompt, and late correlators for each channel
        * will get two outputs from this step
