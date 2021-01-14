# Basic Usage

## Defining Signal Types

For both signal simulation and processing, the `SignalType` must be defined. Currently, `GNSSTools` can simulate and process L1 C/A and L5 signals by default. See the `Custom Signal Types` section below for instructions on defining custom signals in `GNSSTools`. Defining the L1 C/A signal type can be done using the following:

```julia
signal_type = define_l1ca_code_type()
```

This defines a `SignalType` struct that contains the codes for the L1 C/A signal and its carrier frequency. This version does not contain databits by default. Whether a `SignalType` contains databits or not, will not affect data processing. However, dual channel processing is not supported. Either the I or Q channel can be processed, but not both. Defining the L1 C/A signal type this way, places all the codes on the I channel. The Q channel is empty. To define an L1 C/A signal type that has databits, the length of time must be specified.

```julia
signal_type = define_l1ca_code_type(1)
```

Other arguments (`sig_freq`, `prns`, and `databits`) can be provided as well. `sig_freq` is the carrier frequency in Hz, while `prns` is a vector of PRNs that should be included in the signal. This does not mean that they will all be simulated. Multi-PRN simulation is still not available in `GNSSTools` yet. Instead, it includes them in the signal type so they can be generated quickly. All available PRNs for the signal are included by default. A dictionary containing the databits for each PRN can be provided to the argument `databits`, if a custom stream of databits is required. If it is not provided, a random databit stream is assigned to each PRN included in the signal type. The process is identical for defining the L5 signal type, which without databits is simply:

```julia
signal_type = define_l5_code_type()
```

It has one extra argument, `channel`, with available arguments `"both"`, `"I"`, and `"Q"`. `"both"` includes both the I and Q channel codes, while `"I"` and `"Q"` include their respective channel codes.

## Simulating Signals

Signal simulations starts with defining a signal, in the form of a `ReplicaSignal` type. Use:

```julia
signal = definesignal(signal_type, f_s, t_length; [optional args])
```

where `f_s` is the sampling frequency in `Hz` and `t_length` is the length of the signal in `seconds`. Optional arguments such as `prn`, `CN0`, etc can be provided as well. See `definesignal` for defaults. Signal parameters can be redefined after they have been already created using:

```julia
signal = definesignal!(signal; prn=new_prn, phi=new_phi, CN0=new_CN0, ... )
```

Note that a signal's sampling frequency, length, and signal type cannot be changed once it has been defined. All other parameters can be changed. Now, the signal has only been setup for simulation. No actual signal exists yet. Use the following to simulate the signal:

```julia
generatesignal!(signal)
```

The simulated data will be stored in `signal.data`. There is also an accompanying time vector, `signal.t`. Only constant or linear Doppler is simulated using this method. To simulate non-linear Doppler, provide a Doppler curve and associated time vector.

```julia
generatesignal!(signal; doppler_curve=doppler_curve, doppler_t=doppler_t)
```

Doing this will cause the values in `data.f_d` and `data.fd_rate` to be ignore. `doppler_curve` is interpolated and used to generate a signal that contains non-linear Doppler.  

## Loading GNSS Data Files

Data can be loaded into `GNSSTools` for processing. `GNSSTools` assumes that the data is raw IQ data with no header information in the data file. Therefore, specific properties of the data need to be provided in order to read the data correctly. Reading data can be done using:

```julia
data = loaddata(data_type, file_name, f_s, f_center, f_gnss, t_length; ... )
```

`data_type` specifies the bit depth of the complex data. The following can be used:

- `sc4`: 4-bit complex
- `sc8`: 8-bit complex
- `sc16`: 16-bit complex
- `sc32`: 32-bit complex
- `sc64`: 64-bit complex

`f_s` is the sampling rate of the data in Hz, while `f_center` and `f_gnss` are the center frequencies of the receiver and signal, respectively, in Hz. `t_length` is the length of the data to load in seconds. The `skip_to` argument can passes to specify the time in the data to start loading from. By default, the beginning of the data loaded is from 0.001 seconds in the file.

## Processing Signals

All `GNSSSignal` types can be processed. These include `ReplicaSignal` and `GNSSData` types. Do so, the `process` function can be used to perform acquisition and tracking for a single PRN. Note that only individual channels can be processed at a time. Dual channel processing is currently not supported.

```julia
results = process(signal, signal_type, prn, channel="I")
```

`signal_type` is the signal type defined in the above sections. Make sure that it reflects the signal you are trying to process. Note that `channel` is optional, but it defaults to `"I"`. See `process` in the API docs to see the other optional arguments.
