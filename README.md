<img width="300" alt="Krang" src="https://github.com/user-attachments/assets/bd42b144-c8b4-45aa-a331-944cfa9e9794">

# Kerr Ray-tracer for Analytic Null Geodesics (KRANG)

This Julia language package that for a accurate and differentiable null geodesics in the Kerr spacetime.

The package is intended mainly for scientic usage for astrophysical observations, and thus, have constrained the observer to lie at asymptotic infinity.
These algorithms mainly follow the formalism of [Gralla & Lupsasca](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.101.044032), with the exception for a new sub image indexing scheme and a regularized time integral definition.
The ray tracing scheme has been optimized for GPU compatibility and automatic differentiability with [Enzyme.jl](https://enzyme.mit.edu/julia/stable/). 
These considerations allow our algorithms to be easily used in Machine Learning applications with the existing julia infrastructure.

## Documentation
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dominic-chang.github.io/Krang.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dominic-chang.github.io/Krang.jl/dev/)
## Repo Status
[![Build Status](https://github.com/dominic-chang/Krang.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dchang10/Krang.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Coverage](https://codecov.io/gh/dominic-chang/Krang.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/dchang10/Krang.jl)
[![status](https://joss.theoj.org/papers/378df5c54cd21e293b92ac692c21c0ed/status.svg)](https://joss.theoj.org/papers/378df5c54cd21e293b92ac692c21c0ed)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13936258.svg)](https://doi.org/10.5281/zenodo.13936258)
# Installation
Launch your Julia session with then type `]` to move into Pkg mode. Once in pkg mode type
```julia
add Krang
```
