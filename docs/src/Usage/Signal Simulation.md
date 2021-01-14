## Signal Simulation

Signals can be simulated by defining a signal first, specifying characteristics such as Doppler frequency, PRN, C/N₀, etc. By default, thermal and phase noise and ADC quantization are included in signals unless otherwise specified. Defining signals can be done doing:

```julia
signal = definesignal(signal_type, f_s, t_length)
```

where `signal_type` can be defined for L1 C/A and L5 signals using:
