# Overview

GNSSTools is a software suite designed for simulating navigation signals. It is able to perform raw IQ signal generation based off preset orbit information for existing and user defined satellite orbit and constellations. Signal processing can be performed on both simulated and real data, in which acquisition and tracking are performed. Signal tracking currently uses a 1ˢᵗ
order PIF for the delay lock loop (DLL) and a 2 or 3-state Kalman filter for the phase lock loop (PLL). 

GNSSTools is capable of simulating and processing the following GNSS signals:

- L1 C/A
- L5 I/Q
- L1C `IN PROGRESS`

Below is a list of signal characteristics that are simulated

- ADC quantization
- Thermal noise
- Phase noise
- **Random Databits**
    + User can supply custom databit stream if they want to
- Linear and non-linear Doppler frequency effects




