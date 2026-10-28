# Wireless Communication Channel Simulation

**Sanaz Mazrae** · Sharif University of Technology · Communication Systems

MATLAB implementation of six wireless channel models to study the effects of multipath propagation, time-varying fading, Doppler spread, and additive Gaussian noise on complex-baseband signals.

## Channel Models

| Channel | Model |
|---|---|
| `Ray_TI_FF` | Time-invariant Rayleigh, flat fading |
| `Ray_TI_FS` | Time-invariant Rayleigh, frequency-selective fading |
| `Ray_TV_FF` | Time-varying Rayleigh, flat fading |
| `Ray_TV_FS` | Time-varying Rayleigh, frequency-selective fading |
| `Ric_TV_FS` | Time-varying Rician, frequency-selective fading |
| `Awgn` | Additive white Gaussian noise |

The frequency-selective models use a three-tap channel impulse response. Time-varying channels use a sum-of-sinusoids fading model; the Rician channel includes a line-of-sight component. All six models incorporate additive complex Gaussian noise.

## MATLAB Implementation

- [`simulate_channel.m`](matlab/simulate_channel.m): Simulates the selected channel using the interface `simulate_channel(ch_input, ch_type, snr, Tm, fd)`; an optional second output returns channel taps and noise information.
- [`run_channel_analysis.m`](matlab/run_channel_analysis.m): Evaluates impulse responses, temporal fading, frequency selectivity, Doppler correlation, and noise effects using QPSK input signals.

To run the analysis, open MATLAB, navigate to the `matlab` directory, and execute:

```matlab
run_channel_analysis
```

The script saves its figures and numerical summaries to `figures/` and `results/` in the repository's local working directory.

## Report

The [project report](Report.pdf) contains the channel models, simulation parameters, plots, numerical results, and discussion.
