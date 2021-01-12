# Basic Usage

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

Signal simulations starts with defining a signal. Use:

```julia
signal = definesignal(signal_type, f_s, t_length)
```


## Custom Signal Types

For both signal simulation and processing, the `SignalType` must be defined. `SignalType` contains `CodeType` structs, which contain information for a given channel. `SignalType` describes the whole signal, such as codes on both the I and Q channels, bandwidth, and carrier frequency. `CodeType` defines the codes for a specific channel (either I or Q). The user would start by defining a `CodeType`.
