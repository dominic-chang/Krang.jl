```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "KRANG.jl"
  text: "Kerr Raytracer for Analytic Null Geodesics"
  tagline: A performant accurate raytracer for problems in the Kerr spacetime.
  image:
    src: /Krang_logo.png
    alt: Krang.jl
  actions:
    - theme: brand
      text: Getting Started
      link: /getting_started
    - theme: alt
      icon: github
      text: View on Github
      link: https://github.com/dominic-chang/Krang.jl
      text: API 
      link: /api

features:
  - icon: <img width="64" height="64" src="https://metal.juliagpu.org/stable/assets/logo.png" />
    title: GPU Compatible
    details: Type stable and type preserving. GPU compatible with CUDA.jl and Metal.jl.
  - icon: <img width="64" height="64" src="https://enzyme.mit.edu/julia/stable/assets/logo.svg" alt="markdown"/>
    title: Automatic Differentiable, Backwards Propagation compatible
    details: Supports automatie differentiation with Enzyme.jl
---
```