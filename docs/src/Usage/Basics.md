# Basic Usage

For both signal simulation and processing, the `SignalType` must be defined. Currently, `GNSSTools` can simulate and process L1 C/A and L5 signals by default. See the `Custom Signal Types` section below for instructions on defining custom signals in `GNSSTools`. Defining the L1 C/A signal type can be done using the following:

```julia
signal_type = define_l1ca_code_type()
```

This defines a `SignalType` struct that contains the codes for the L1 C/A signal and its carrier frequency. This version does not contain databits by default. Whether a `SignalType` contains databits or not, will not affect data processing. However, dual channel processing is not supported. Either the I or Q channel can be processed, but not both. To define an L1 C/A signal type that has databits, the length of time must be specified.

```julia
signal_type = define_l1ca_code_type(1)
```

##

## Custom Signal Types

For both signal simulation and processing, the `SignalType` must be defined. `SignalType` contains `CodeType` structs, which contain information for a given channel. `SignalType` describes the whole signal, such as codes on both the I and Q channels, bandwidth, and carrier frequency. `CodeType` defines the codes for a specific channel (either I or Q). The user would start by defining a `CodeType`.
