<img width="300" alt="Krang" src="https://github.com/user-attachments/assets/bd42b144-c8b4-45aa-a331-944cfa9e9794">

# Kerr Ray-tracer for Analytic Null Geodesics (KRANG)

This Julia package provides accurate, differentiable null geodesics in the Kerr spacetime.

The package is intended mainly for scientific use in astrophysical observations, so the observer is constrained to lie at asymptotic infinity.
These algorithms mainly follow the formalism of [Gralla & Lupsasca](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.101.044032), with the exception of a new sub-image indexing scheme and a regularized time integral definition.
The ray tracing scheme has been optimized for GPU compatibility and automatic differentiability with [Enzyme.jl](https://enzyme.mit.edu/julia/stable/). 
These considerations allow the algorithms to be used in machine-learning applications with the existing Julia infrastructure.

## Documentation
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dominic-chang.github.io/Krang.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dominic-chang.github.io/Krang.jl/dev/)

## Repo Status
[![Build Status](https://github.com/dominic-chang/Krang.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dominic-chang/Krang.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Coverage](https://codecov.io/gh/dominic-chang/Krang.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/dominic-chang/Krang.jl)
[![status](https://joss.theoj.org/papers/378df5c54cd21e293b92ac692c21c0ed/status.svg)](https://joss.theoj.org/papers/378df5c54cd21e293b92ac692c21c0ed)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13936258.svg)](https://doi.org/10.5281/zenodo.13936258)

# Installation
Launch a Julia session, then type `]` to enter Pkg mode. Once in Pkg mode, type:
```julia
add Krang
```

# Simple Benchmark

Ray tracing the coordinates of equatorial emitters to pixels on the observer's screen. (Apple M1 Pro, 8 performance cores)

| Param. | Value |
|---|---|
| Spin | a = 0.94 |
| Observer Inclination | θ₀ = 45° |
| Screen Size | ±15 M |
| Disc Extent | 2.024 – 20 M |
| Observer Distance | r = 1000 M |
| Sub-Image | n = 0 direct |

## CPU time

Single thread execution time returning emission
radius, azimuth and time.

| Code | Family | Per pixel |
|---|---|---:|
| `Krang` | analytic | 2.86 µs |
| `kgeo` | analytic | 2.95 µs |
| `aart` | analytic | 5.01 µs |
| `ipole` | numerical | 2063 µs |

## Clock time

**note**: `kgeo` and `aart` are numpy-vectorized with no multithreaded entry.

| Code | Family | 64² | 128² | 256² | 256² / 8 threads | Speedup |
|---|---|---:|---:|---:|---:|---:|
| `Krang` | analytic | 11.5 ms | 46.6 ms | 187.3 ms | 29.0 ms | 6.5× |
| `kgeo` | analytic | 13.2 ms | 48.1 ms | 193.3 ms | — | — |
| `aart` | analytic | 21.8 ms | 83.1 ms | 328.4 ms | — | — |
| `ipole` | numerical | 8.31 s | 32.83 s | 135.21 s | 21.08 s | 6.4× |

