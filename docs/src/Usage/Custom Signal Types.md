## Custom Signal Types

For both signal simulation and processing, the `SignalType` must be defined. `SignalType` contains `CodeType` structs, which contain information for a given channel. `SignalType` describes the whole signal, such as codes on both the I and Q channels, bandwidth, and carrier frequency. `CodeType` defines the codes for a specific channel (either I or Q). The user would start by defining a `CodeType`.

```julia

```
