---
title: 'COSMOS: A numerical relativity code specialized for PBH formation'
tags:
  - C++
  - cosmology
  - dynamics
  - black hole formation
  - gravitational collapse
authors:
  - name: Chul-Moon Yoo
    orcid: 0000-0002-9928-4757
    # equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Albert Escrivà
    # equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Tomohiro Harada
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - name: Hayami Iizuka
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - name: Taishi Ikeda
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 4
  - name: Yasutaka Koga
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 5
  - name: Hirotada Okawa
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 6
  - name: Daiki Saito
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Masaaki Shimada
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Koichiro Uehara
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
#   - given-names: Ludwig
    # dropping-particle: van
    # surname: Beethoven
    # affiliation: 3
affiliations:
 - name: Graduate School of Science, Nagoya University, Japan
   index: 1
#    ror: 00hx57361
 - name: National Astronomical Observatory of Japan (NAOJ), Japan
   index: 2
 - name: Department of Physics, Rikkyo University, Japan
   index: 3
 - name: Niels Bohr International Academy, Niels Bohr Institute, Denmark
   index: 4
 - name: Yukawa Institute of Theoretical Physics, Kyoto University, Japan
   index: 5
 - name: Faculty of Software and Information Technology, Aomori University, Japan
   index: 6
date: 5 February 2025
bibliography: paper.bib


# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---


# Summary


Primordial black holes (PBHs) are those generated in the early universe without experience of the form of a star.
It has been pointed out that PBHs may be candidates for black holes and compact objects of various masses in the universe
or a major component of dark matter, and PBHs are attracting much attention. In the standard formation process, PBHs are formed from super-horizon primordial fluctuation with a non-linearly large initial amplitude.
In order to follow the whole non-linear gravitational dynamics, one has to rely on numerical relativity solving Einstein equations numerically.
`COSMOS` [@Yoo:2013yea; @Okawa:2014nda] and `COSMOS-S` [@Yoo:2021fxs] provide simple tools for the simulation of PBH formation.
`COSMOS` and `COSMOS-S` are C++ packages for solving Einstein equations in 3+1 dimension and spherical symmetry (1+1 dimension), respectively.
It was originally translated from SACRA code [@Yamamoto:2008js] into C++ and developed specialized for PBH formation.  
<!-- The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration). -->


# Statement of need


In the simulation of PBH formation, since there is a hierarchy between the size of the collapsing region and cosmological expansion scale, an efficient resolution refinement procedure is needed. In order to resolve the collapsing region, non-Cartesian scale-up coordinates [@Yoo:2018pda] and a fixed mesh-refinement procedure are implemented. To achieve a practically acceptable computational speed, an OpenMP package is used for the parallelization. No other packages are required, and the functionality is minimal. Therefore it would be easy to use for beginners of numerical relativity. Perfect fluid with linear equation of states and massless scalar field are implemented as matter fields. Once users understand the source code to some extent, the system can be easily extended to various scientifically interesting settings.




<!--
`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).


`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike. -->


# Physical settings


Einstein equations
<!-- \begin{equation} -->
$$
G_{\mu\nu}=R_{\mu\nu}-\frac{1}{2}Rg_{\mu\nu}=\frac{8\pi G}{c^4}T_{\mu\nu}
$$
<!-- \end{equation} -->
are solved, where $G_{\mu\nu}$, $g_{\mu\nu}$, $R_{\mu\nu}$, $R$, $G$, $c$ and $T_{\mu\nu}$
are the Einstein tensor, metric tensor, Ricci tensor, Ricci scalar, Newtonian gravitational constant, speed of light and energy momentum tensor, respectively.
The energy momentum tensor can be divided into the fluid and scalar field contributions as follows:
<!-- \begin{equation} -->
$$
T_{\mu\nu}=T^{\rm SC}_{\mu\nu}+T^{\rm FL}_{\mu\nu}
$$
<!-- \end{equation} -->
with
<!-- \begin{equation} -->
$$
T^{\rm SC}_{\mu\nu}=\nabla_\mu\phi\nabla_\nu\phi-\frac{1}{2}g_{\mu\nu}\nabla^\lambda\phi\nabla_\lambda\phi
$$
<!-- \end{equation} -->
and
<!-- \begin{equation} -->
$$
T^{\rm FL}_{\mu\nu}=(\rho+P)u_\mu u_\nu+Pg_{\mu\nu},
$$
<!-- \end{equation} -->
where $\nabla$, $\phi$, $\rho$, $u_\mu$ and $P$ are the covariant derivative for $g_{\mu\nu}$, scalar field, fluid energy density and pressure, respectively.
The equations of motion for the scalar field
<!-- \begin{equation} -->
$$
\nabla^\mu\nabla_\mu \phi=0
$$
<!-- \end{equation} -->
and the fluid
<!-- \begin{equation} -->
$$
\nabla^\mu T^{\rm FL}_{\mu\nu}=0
$$
<!-- \end{equation} -->
are also solved.
Readers are asked to refer to standard textbooks of numerical relativity (e.g., @gourgoulhon20123+1; @shibata2016numerical) for how to rewrite these equations into the form suitable for numerical integration.


As for the initial data, we adopt the long-wavelength growing-mode solutions
up through the next-leading order of the expansion parameter $\epsilon=k(aH_b)\ll1$,
where $1/k$ gives the characteristic comoving scale of the inhomogeneity, and $a$ and $H_b$ are the scale factor and Hubble expansion rate in the reference background universe.
The initial data can be characterized by a function of the spatial coordinates $\vec x$ as the curvature perturbation $\zeta(\vec x)$ for adiabatic fluctuations [@Harada:2015yda; @Yoo:2024lhp; @Yoo:2020lmg] and iso-curvature perturbation $\Upsilon(\vec x)$ for
massless scalar iso-curvature [@Yoo:2021fxs]. Since the space is filled with fluid, the initial fluid distribution can be generated by imposing the constraint equations.
Then the Hamiltonian and momentum constraints are initially satisfied within the machine's precision. Therefore we do not solve the constraint equations as in the case of asymptotically flat systems with the existence of a vacuum region. Elliptic solvers for constraint equations are not included in this package for the above reasons.


# Examples
<!--
Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format. -->
Several examples are included in the package.
These are just for demonstration and exercise for users, and we do not care about the precision of the examples. The resolution is kept to a minimum. In the figures below, the length scale is normalized by the size $L$ of the box for the numerical simulation.


### COSMOS (3+1 dimensional simulation)


- Evolution of a single mode perturbation


The evolution of sinusoidal small fluctuation is given as an example, which can be compared with the corresponding linear perturbation (see \autoref{fig:kap}).


![The time evolution of the trace of the extrinsic curvature tr$K$ is compared with the solution of the linear perturbation equation. $H$ in the vertical axis label is the Hubble expansion rate in the background universe model.\label{fig:kap}](kap.pdf){height="6cm"}




- Adiabatic spherically symmetric initial fluctuation


The scalar field is absent in this example. The setting is similar to that in @Yoo:2020lmg.
We also attach the data obtained by solving the Einstein equations until an apparent horizon is found (see \autoref{fig:alp} and \autoref{fig:AH}).


![The lapse function (``$tt$-component" of the metric) on the $x$-axis at the time when an apparent horizon is found.\label{fig:alp}](alp_2d.pdf){height="8cm"}


![The shape of the apparent horizon when it is found.\label{fig:AH}](AH_tex.pdf){height="6cm"}


- Spherically symmetric iso-curvature


The setting is similar to that in @Yoo:2021fxs.
We also attach the data obtained by solving the Einstein equations until an apparent horizon is found.






### COSMOS-S (spherically symmetric simulation)


- Adiabatic spherically symmetric initial fluctuation


The physical parameter setting is the same as the corresponding example for 3+1 dimensional simulation. But resolution is finer in this example of the spherically symmetric 1+1 code.


- Spherically symmetric iso-curvature


The physical parameter setting is the same as the corresponding example for 3+1 dimensional simulation. But resolution is finer in this example of the spherically symmetric 1+1 code.


- Type II-B PBH formation


PBH formation from adiabatic fluctuation with an extremely large initial amplitude is given as an example. The setting is similar to that in @Uehara:2024yyp. One can find the non-trivial trapping horizon configuration as \autoref{fig:horizon}.


![Trapping horizon trajectories.\label{fig:horizon}](horizon.pdf){height="6cm"}


<!--
# Figures


Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.


Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% } -->


# Acknowledgements


A.E. acknowledges support from the JSPS Postdoctoral Fellowships for Research in Japan (Graduate School of Sciences, Nagoya University).
D.S. is supported in part by JSPS KAKENHI Grant No. 24KJ1223. K.U. would like to
take this opportunity to thank the “THERS Make New Standards Program for the Next
Generation Researchers” supported by JST SPRING, Grant Number JPMJSP2125. C.Y. is
supported in part by JSPS KAKENHI Grant Nos. 20H05850, 20H05853 and 24K07027.




# References



