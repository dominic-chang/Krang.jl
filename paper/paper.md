---
title: 'Krang: Kerr Raytracer for Analyitc Null Geodesics'
tags:
  - julia
  - black hole 
  - general relativity
  - astronomy
  - raytracing
authors:
  - name: Dominic Chang
    orcid: 0000-0001-9939-5257
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Physics, Harvard University
   index: 1
 - name: Black Hole Initiative at Harvard University
   index: 2
date: 05 Sep 2024
bibliography: paper.bib

---

# Summary
`Krang` is a Julia[^1] [@Bezanson2015] package that implements efficient algorithms for raytracing emission geometries in the Kerr black hole space time.
It is GPU compatible and is specialized for studies of sub-image contributions from gravitational lensed sources [@JohnsonRing]. 
Such algorithms are of interest for modeling the sources seen by Very Long Baseline Interferometry (VLBI) observations of Low Luminosity Active Galactic Nuclei (LLAGN) such as those imaged by the Event Horizon Telescope Collaboration (EHTC).

[^1]: https://julialang.org

# Statement of need

Studies of the near horizon scale science around supermassive black holes have increased in interest due to results from the horizon scale observations of supermassive black holes by the Gravity Collaboration, the Atacama Large Millimeter Array (ALMA), and the EHTC.
The Event Horizon Telescope, in particular [EHT, @M87PaperII], produced the first images of the shadows of the supermassive black holes in the centres of M87 [@M87PaperI] and the Milky Way [@SgrAPaperI] galaxies at event horzon scales.
The emission from these sources originating from relativistic plasmas that are accreted by supermassive black holes in the presence of a magnetic field.
Scientific analysis of the data from these sources often require complicated source modelling that includes various relativistic effects which leave characteristic signatures in the observed images.
The large scale of the black holes allow for many of these effects to be described within the geometric optics limit of electro-magnetism. Ray tracing techniques thus present viable options for modeling images of super massive black holes.

A relativistic image feature that has been theorized to exist, but is yet to be observed, are the individual sub-image contributions to the overall image structure known as photon-rings, [@JohnsonRing].
Photon rings are of particular interest because of their strong dependence on gravitational effects and their insensitivity to variations in the emission physics around the black hole.
The observation of a photon ring would therefore serve as effective probe for measurements of black hole characteristics like spin, or help facilitate tests of gravity.
This feature could potentially be seen in the near future with a recent as a space extension to the EHT aimed of detection and measurement [@lupsasca2024blackholeexplorerphoton].

Scientific studies of black hole images within the current software landscape often requires great compromise due to the computational complexity of the problem, and the incompatibility of the existing software with other analysis packages.
It is difficult, for example, to apply machine learning optimization algorithms to existing raytracing codes in python implementations since they typically rely on special functions that are not implemented within existing machine learning frameworks.
Existing python implementations are also bounded to CPU evaluations, limiting their capability of accessing acceleration from specialized hardware.
A Julia implementation of raytracing algorithms is thus beneficial because of the languages modular design, differentiable programming and efficient execution.
Krang therefore benefits from the 'plug and play' nature of the julia programming language, allowing for easy development and synergy with existing analysis pipelines, and fast CPU/GPU executions.


# Similar Packages
- `AART` [@cardenas-avendano]: An Adaptive Analytical Ray Tracing code for geodesics in the Kerr space time in python.

- `KerrBam` [@Palumbo_2022]: An adaptive analytical raytracing code for equatorial synchrotron models in python.

- `Gradus` [@Gradus]: A Julia implementation of a numeric, differentiable, general relativistic raytracer.

# Acknowledgements
The authors thank Iniyan Natarajan, Alexander Plavin, and Paul Tiede for helpful discussions. Support for this work was provided by the NSF (AST-1935980, AST-2034306) and by the Gordon and Betty Moore Foundation through grant GBMF-10423. This work was supported by the Black Hole Initiative, which is funded by grants from the John Templeton Foundation (Grant #62286) and the Gordon and Betty Moore Foundation (Grant GBMF-8273), although the opinions expressed in this work are those of the author(s) and do not necessarily reflect the views of these Foundations.

# References
