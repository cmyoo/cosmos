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
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Albert Escrivà
    orcid: 0000-0001-5483-8034
    # equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 3
  - name: Tomohiro Harada
    orcid: 0000-0002-9085-9905
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 4
  - name: Hayami Iizuka
    orcid: 0009-0001-6604-7763
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 4
  - name: Taishi Ikeda
    orcid: 000-0002-9076-1027
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 5
  - name: Yasutaka Koga
    orcid: 0000-0002-9579-5787
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 6
  - name: Hirotada Okawa
    orcid: 0000-0001-7372-5131
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 7
  - name: Daiki Saito
    orcid: 0000-0003-1624-9268
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Masaaki Shimada
    orcid: 0009-0001-2144-575X
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Koichiro Uehara
    orcid: 0009-0006-3039-6829
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
#   - given-names: Ludwig
    # dropping-particle: van
    # surname: Beethoven
    # affiliation: 3
affiliations:
 - name: Graduate School of Science, Nagoya University, Japan
   index: 1
 - name: Kobayashi Maskawa Institute, Nagoya University, Japan
   index: 2
#    ror: 00hx57361
 - name: National Astronomical Observatory of Japan (NAOJ), Japan
   index: 3
 - name: Department of Physics, Rikkyo University, Japan
   index: 4
 - name: Center of Gravity, Niels Bohr Institute, Denmark
   index: 5
 - name: Yukawa Institute for Theoretical Physics, Kyoto University, Japan
   index: 6
 - name: Faculty of Software and Information Technology, Aomori University, Japan
   index: 7
date: 5 February 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Primordial black holes (PBHs) are black holes generated in the early universe without experience of the form of a star.
It has been pointed out that PBHs may be candidates for black holes and compact objects of various masses in the universe
or a major component of dark matter. 
In particular, PBHs have been attracting much attention in the recent development of gravitational wave observation. 
In the standard formation process, PBHs are formed from super-horizon primordial fluctuations with non-linearly large initial amplitude.
In order to follow the whole non-linear gravitational dynamics, one has to rely on numerical relativity solving Einstein equations.
`COSMOS` [@Yoo:2013yea; @Okawa:2014nda] and `COSMOS-S` [@Yoo:2021fxs] provide simple tools for the simulation of PBH formation. 
`COSMOS` and `COSMOS-S` are C++ packages for solving Einstein equations in 3+1 dimensions and spherical symmetry (1+1 dimensions), respectively. 
It was originally translated from SACRA code [@Yamamoto:2008js] into C++ and has been developed specialized for PBH formation. 
In this paper, we do not describe all scientific results obtained by using `COSMOS` or `COSMOS-S`. 
The readers who are interested in the past achievments may refer to @Yoo:2013yea; @Okawa:2014nda; @Yoo:2014boa; @Okawa:2014sxa; @Ikeda:2015hqa; @Brito:2015yga; @Brito:2015yfh; @Okawa_2015; @Yoo:2016kzu; @Yoo:2018pda; @Yoo:2021fxs; @Yoo:2024lhp; @Escriva:2024lmm; @Escriva:2024aeo; @Shimada:2024eec [^1]. 

[^1]: In most of the references: @Yoo:2013yea; @Okawa:2014nda; @Yoo:2014boa; @Okawa:2014sxa; @Ikeda:2015hqa; @Brito:2015yga; @Brito:2015yfh; @Okawa_2015; @Yoo:2016kzu; @Yoo:2018pda; @Yoo:2021fxs; @Yoo:2024lhp; @Escriva:2024lmm; @Escriva:2024aeo; @Shimada:2024eec, 
additional functions and packages have been implemented to meet the requirements for individual settings. 
Therefore the results may not be obtained by simply running the public code. 

# Statement of need

In the simulation of PBH formation, since there is a hierarchy between the size of the collapsing region and cosmological expansion scale, an efficient resolution refinement procedure is needed. 
In order to resolve the collapsing region, non-Cartesian scale-up coordinates [@Yoo:2018pda] and a fixed mesh-refinement procedure [@Yoo:2024lhp] are implemented in `COSMOS`. 
The 1+1 dimensional simulation code `COSMOS-S` [@Yoo:2021fxs] is derived from `COSMOS` with the CARTOON method [@Alcubierre:1999ab]. 
To achieve a practically acceptable computational speed, an OpenMP package is used for the parallelization. 
No other packages are required, and the functionality is minimal. 
Therefore it would be easy to use for beginners of numerical relativity. 
A perfect fluid with a linear equation of states and a massless scalar field are implemented as matter fields. 
Once users understand the source code to some extent, the system can be easily extended to various scientifically interesting settings.


# Physical system settings

Einstein equations
<!-- \begin{equation} -->
$$
G_{\mu\nu}=R_{\mu\nu}-\frac{1}{2}Rg_{\mu\nu}=\frac{8\pi G}{c^4}T_{\mu\nu}
$$
<!-- \end{equation} -->
are solved, where $G_{\mu\nu}$, $g_{\mu\nu}$, $R_{\mu\nu}$, $R$, $G$, $c$ and $T_{\mu\nu}$
are the Einstein tensor, metric tensor, Ricci tensor, Ricci scalar, Newtonian gravitational constant, speed of light and energy momentum tensor, respectively.
The energy momentum tensor is divided into the fluid and scalar field contributions as follows:
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
where $\nabla$, $\phi$, $\rho$, $u_\mu$ and $P$ are the covariant derivative associated with $g_{\mu\nu}$, scalar field, fluid energy density, four-velocity and pressure, respectively. 
The pressure and the energy density are assumed to satisfy the linear equation of state $P=w\rho$ with $w$ being a constant. 
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
Readers are asked to refer to standard textbooks of numerical relativity (e.g., @gourgoulhon20123+1; @shibata2016numerical) to learn how to rewrite these equations into a form suitable for numerical integration. 
To solve the fluid equations of motion, we basically follow the scheme discussed in @Kurganov:2000ovy; @Shibata:2005jv. 

As for the initial data, we adopt the long-wavelength growing-mode solutions
up through the next-leading order of the expansion parameter $\epsilon=k/(aH)\ll1$,
where $1/k$ gives the characteristic comoving scale of the inhomogeneity, and $a$ and $H$ are the scale factor and Hubble expansion rate in the reference background universe.
The initial data can be characterized by a function of the spatial coordinates $\vec x$ as the curvature perturbation $\zeta(\vec x)$ for adiabatic fluctuations [@Harada:2015yda; @Yoo:2024lhp; @Yoo:2020lmg] and iso-curvature perturbation $\Upsilon(\vec x)$ for
massless scalar iso-curvature [@Yoo:2021fxs]. 
Since the space is filled with the fluid, the initial fluid distribution can be generated to meet the constraint equations included in the Einstein equations. 
Then, the constraint equations are initially satisfied within the machine's precision. 
Therefore, the constraint equations are not solved by integrating elliptic differential equations. 
This is very different from the standard method to obtain the initial data for spacetimes 
with asymptotically flat vacuum regions. 
This is why elliptic solvers for constraint equations are not included in this package.

# Examples

Several examples are included in the package.
These examples are intended primarily for demonstration and instructional purposes, and thus, the precision is not a primary concern. The resolution has been intentionally kept to a minimum.
In the figures below, the length scale is normalized by the size $L$ of the box for the numerical simulation. 


### COSMOS (3+1 dimensional simulation)

- Evolution of a single-mode perturbation

The evolution of sinusoidal small fluctuation is given as an example, which can be compared with the corresponding linear perturbation (see \autoref{fig:kap}).

![The time evolution of the trace of the extrinsic curvature tr$K$ is compared with the solution of the linear perturbation equation.\label{fig:kap}](kap.pdf){height="6cm"}



- Adiabatic spherically symmetric initial fluctuation

The scalar field is absent in this example. The setting is similar to that in @Yoo:2020lmg.
We also attach the data obtained by solving the Einstein equations until an apparent horizon is found (see \autoref{fig:alp} and \autoref{fig:AH}).

![The lapse function (``$tt$-component" of the metric) on the $xy$-plane at the time when an apparent horizon is found.\label{fig:alp}](alp_2d.pdf){height="8cm"}

![The shape of the apparent horizon when it is found.\label{fig:AH}](AH_tex.pdf){height="6cm"}

- Spherically symmetric iso-curvature

The setting is similar to that in @Yoo:2021fxs.
We also attach the data obtained by solving the Einstein equations until an apparent horizon is found.





### COSMOS-S (spherically symmetric simulation)

- Adiabatic spherically symmetric initial fluctuation

The physical parameter setting is the same as the corresponding example for the 3+1 dimensional simulation. However, 
the resolution is finer in this example of the spherically symmetric 1+1 code.

- Spherically symmetric iso-curvature

The physical parameter setting is the same as the corresponding example for the 3+1 dimensional simulation. However, 
the resolution is finer in this example of the spherically symmetric 1+1 code.

- Type II-B PBH formation

PBH formation from adiabatic fluctuation with extremely large initial amplitude is given as an example. The setting is similar to that in @Uehara:2024yyp. One can find the non-trivial trapping horizon configuration as \autoref{fig:horizon}.

![Trapping horizon trajectories.\label{fig:horizon}](horizon.pdf){height="6cm"}


# Acknowledgements

A.E. acknowledges support from the JSPS Postdoctoral Fellowships for Research in Japan (Graduate School of Sciences, Nagoya University).
K.U. would like to take this opportunity to thank the “THERS Make New Standards Program for the Next
Generation Researchers” supported by JST SPRING, Grant Number JPMJSP2125. 
T.I. acknowledges support from VILLUM Foundation (grant no. VIL37766) and the DNRF Chair program (grant no. DNRF162) by the Danish National
Research Foundation. 
This work is supported in part by JSPS KAKENHI Grant Nos. 20H05850 (C.Y.), 20H05853 (T.H., C.Y.), 21K20367 (Y.K.), 23KK0048 (Y.K.), 24K07027 (T.H., C.Y.) and 24KJ1223 (D.S.).


# References


