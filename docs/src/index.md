```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "KRANG.jl"
  text: "Kerr Ray-tracer for Analytic Null Geodesics"
  tagline: A differentiable and accurate ray tracer for problems in the Kerr spacetime.
  image:
    src: /Krang_logo.png
    alt: Krang.jl
  actions:
    - theme: brand
      text: Getting Started
      link: /what_is_krang
    - theme: alt
      text: View on Github
      link: https://github.com/dominic-chang/Krang.jl
    - theme: alt
      text: API 
      link: /api

features:
  - icon: <img width="64" height="64" src="https://metal.juliagpu.org/stable/assets/logo.png" />
    title: GPU Compatible
    details: Type stable and type preserving. GPU compatible with CUDA.jl and Metal.jl.
  - icon: <img width="64" height="64" src="https://enzyme.mit.edu/julia/stable/assets/logo.svg" alt="markdown"/>
    title: Automatic Differentiable, Backwards Propagation compatible
    details: Supports automatic differentiation with Enzyme.jl
---
```