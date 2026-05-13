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
  - name: Hirotada Okawa
    orcid: 0000-0001-7372-5131
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - name: Albert Escrivà
    orcid: 0000-0001-5483-8034
    # equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Tomohiro Harada
    orcid: 0000-0002-9085-9905
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 4
  - name: Hayami Iizuka
    orcid: 0009-0001-6604-7763
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 4
  - name: Taishi Ikeda
    orcid: 0000-0002-9076-1027
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 5
  - name: Yasutaka Koga
    orcid: 0000-0002-9579-5787
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 6
  - name: Daiki Saito
    orcid: 0000-0003-1624-9268
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 7
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
 - name: Kobayashi-Maskawa Institute for the Origin of Particles and the Universe, Nagoya University, Japan
   index: 2
 - name: Faculty of Software and Information Technology, Aomori University, Japan
   index: 3
#    ror: 00hx57361
# - name: National Astronomical Observatory of Japan (NAOJ), Japan
#   index: 3
 - name: Department of Physics, Rikkyo University, Japan
   index: 4
 - name: Center of Gravity, Niels Bohr Institute, Denmark
   index: 5
 - name: Department of Information and Computer Science, Osaka Institute of Technology, Japan
   index: 6
 - name: Department of Physics, Kyoto University, Japan
   index: 7
date: 5 February 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Primordial black holes (PBHs) are black holes generated in the early universe without having gone through stellar evolution.
It has been hypothesized that PBHs may be candidates for black holes and compact objects of various masses in the universe
or a major component of dark matter. 
In particular, PBHs have been attracting much attention in the recent development of gravitational wave observation. 
In the standard formation process, PBHs are formed from super-horizon primordial fluctuations with non-linearly large initial amplitude.
In order to simulate the non-linear gravitational dynamics of PBH formation, one has to rely on numerical relativity solvers to approximate the solution of the Einstein equations.
`COSMOS` [^1] [@Yoo:2013yea; @Okawa:2014nda] provides <!-- and `COSMOS-S` [@Yoo:2021fxs] provide  -->
 simple tools for the simulation of PBH formation (see `COSMOS-S` [^2] for a spherically symmetric version of `COSMOS`, which is not discussed in this paper). 
`COSMOS` is a <!-- and `COSMOS-S` are  -->
C++ package for solving the Einstein equations in 3+1 dimensions. <!-- and spherical symmetry (1+1 dimensions), respectively.  -->
It was originally translated from SACRA code [@Yamamoto:2008js] into C++ and has been developed specifically for the simulation of PBH formation [^3]. 
Past publications that use COSMOS for simulation include @Yoo:2013yea; @Okawa:2014nda; @Yoo:2014boa; @Okawa:2014sxa; @Ikeda:2015hqa; @Brito:2015yga; @Brito:2015yfh; @Okawa_2015; @Yoo:2016kzu; @Yoo:2018pda; @Yoo:2024lhp; @Escriva:2024lmm; @Escriva:2024aeo [^4]. <!-- @Yoo:2021fxs; @Shimada:2024eec  -->

[^1]: https://github.com/cmyoo/cosmos .
[^2]: https://github.com/cmyoo/cosmos-s .
[^3]: C.Y. and H.O are the main contributors of this code, and other authors used the numerical code during the development and operational stages and contributed in part to its development and improvement.
[^4]: In these works, additional functions and packages have been implemented that may not appear in the public release of COSMOS. Therefore the results may not be obtained by simply running the public code. 

# Statement of need

In the simulation of PBH formation, the presence of multiple lengthscales (the size of the collapsing region and that of cosmological expansion) necessitates an efficient resolution refinement procedure. 
In order to resolve the collapsing region, non-Cartesian scale-up coordinates [@Yoo:2018pda] and a fixed mesh-refinement procedure [@Yoo:2024lhp] are implemented in `COSMOS`. <!-- The 1+1 dimensional simulation code `COSMOS-S` [@Yoo:2021fxs] is derived from `COSMOS` with the CARTOON method [@Alcubierre:1999ab].  -->
In its model, COSMOS uses a perfect fluid with a linear equation of state and a massless scalar field as matter fields. 
To achieve a practically acceptable computational speed, OpenMP is used for the parallelization. 
COSMOS has no other dependencies, which makes for an easier installation. 
Once users understand the source code to some extent, the system can be easily extended to various scientifically interesting settings.


# Physical system settings

At the core of COSMOS, the Einstein equations
<!-- \begin{equation} -->
$$
G_{\mu\nu}=R_{\mu\nu}-\frac{1}{2}Rg_{\mu\nu}=\frac{8\pi G}{c^4}T_{\mu\nu}
$$ 
<!-- \end{equation} --> 
are solved, where $G_{\mu\nu}$, $g_{\mu\nu}$, $R_{\mu\nu}$, $R$, $G$, $c$ and $T_{\mu\nu}$
are the Einstein tensor, the metric tensor, the Ricci tensor, the Ricci scalar, the Newtonian gravitational constant, the speed of light, and the energy-momentum tensor, respectively. 
The energy-momentum tensor is divided into fluid and scalar field contributions as
<!-- \begin{equation} -->
$$
T_{\mu\nu}=T^{\rm SC}_{\mu\nu}+T^{\rm FL}_{\mu\nu}, 
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
where $\nabla$ denotes taking the covariant derivative using $g_{\mu\nu}$, 
$\phi$ is the scalar field, and $\rho$, $u_\mu$ and $P$ are the energy density, the four-velocity, and the pressure of the fluid, respectively. 
The pressure and the energy density are assumed to satisfy the linear equation of state $P=w\rho$ with $w$ being a constant. 
The equations of motion for the scalar field
<!-- \begin{equation} -->
$$
\nabla^\mu\nabla_\mu \phi=0, 
$$
<!-- \end{equation} -->
and the fluid
<!-- \begin{equation} -->
$$
\nabla^\mu T^{\rm FL}_{\mu\nu}=0, 
$$
<!-- \end{equation} -->
are also solved.
Readers are asked to refer to standard textbooks of numerical relativity (e.g., @gourgoulhon20123+1; @shibata2016numerical) on how to rewrite these equations into a form suitable for numerical integration. 
To solve the fluid equations of motion, we basically follow the scheme discussed in @Kurganov:2000ovy; @Shibata:2005jv. 

As for the initial data, we adopt the long-wavelength growing-mode solutions
up to (and including) the next-leading order of the expansion parameter $\epsilon=k/(aH)\ll1$,
where $1/k$ gives the characteristic comoving scale of the inhomogeneity, and $a$ and $H$ are the scale factor and Hubble expansion rate in the reference background universe.
The initial data can be characterized by a function of the spatial coordinates $\vec x$ as the curvature perturbation $\zeta(\vec x)$ for adiabatic fluctuations [@Harada:2015yda; @Yoo:2024lhp; @Yoo:2020lmg] and iso-curvature perturbation $\Upsilon(\vec x)$ for
massless scalar iso-curvature [@Yoo:2021fxs]. 
Since the space is filled with the fluid, the initial fluid distribution can be generated to meet the constraint equations included in the Einstein equations. 
Therefore the constraint equations are initially satisfied to within machine precision, and need not be solved by integrating elliptic differential equations. 
This approach differs from the standard method of obtaining the initial data for spacetimes with asymptotically flat vacuum regions, and is the reason why elliptic solvers are not included in COSMOS. 

# Examples

Three examples are included in the package.
These examples are intended primarily for demonstration and instructional purposes, and thus, and thus they are not highly accurate. 
The resolution has been intentionally kept as low as possible. 
In the figures below, the length scale is normalized by the size $L$ of the box for the numerical simulation. 


<!-- ### COSMOS (3+1 dimensional simulation) -->

### Evolution of a single-mode perturbation

The evolution of a small sinusoidal fluctuation is given as an example, which can be compared to the corresponding linear perturbation (see \autoref{fig:kap}).

![The time evolution of the trace of the extrinsic curvature tr$K$ is compared with the solution of the linear perturbation equation.\label{fig:kap}](kap.pdf){height="6cm"}



### Adiabatic spherically symmetric initial fluctuation

Here, we consider the adiabatic perturbation generated by the initial curvature perturbation without scalar field contribution. 
As is described in @Harada:2015yda, once the spatial profile of the curvature perturbation is specified,
the growing mode solution can be described in the long-wavelength approximation. 
More details including the specific functional form of the curvature perturbation can be found in [^5] (see also @Yoo:2020lmg). 
In the repository, we include the data file `ini_all.dat` necessary to reconstruct the geometry and matter distribution at the time an apparent horizon is found (see \autoref{fig:alp} and \autoref{fig:AH} which can be generated by following the instructions in the example code).  

![The lapse function (``$tt$-component" of the metric) on the $xy$-plane at the time when an apparent horizon is found. The blue, green, and purple meshes show the region covered by the lowest, 1st, and 2nd mesh refinement layers, respectively.\label{fig:alp}](alp_2d.pdf){height="7cm"}

![The shape of the apparent horizon when it is found.\label{fig:AH}](AH_tex.pdf){height="5cm"}

[^5]: https://github.com/cmyoo/cosmos/wiki/Adiabatic-spherically-symmetric-initial-fluctuation .

### Spherically symmetric iso-curvature

Here, we consider the iso-curvature perturbation generated by a massless scalar field.
We assume that the massless scalar field does not contribute to the background metric of the gradient expansion.
Then, as is described in @Yoo:2021fxs, once the 
spatial profile of the scalar field is specified at the leading order of the gradient expansion,
the growing mode solution can be described in the long-wavelength approximation. 
More details including the specific functional form of the curvature perturbation can be found in [^6] (see also @Yoo:2021fxs).
In the repository, we include the data file `ini_all.dat` necessary to reconstruct the geometry and matter distribution at the time an apparent horizon is found, as for the adiabatic case. 

[^6]: https://github.com/cmyoo/cosmos/wiki/Spherically-symmetric-isocurvature .


<!--
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

-->

# Acknowledgements

A.E. acknowledges support from the JSPS Postdoctoral Fellowships for Research in Japan (Graduate School of Sciences, Nagoya University).
K.U. would like to take this opportunity to thank the “THERS Make New Standards Program for the Next
Generation Researchers” supported by JST SPRING, Grant Number JPMJSP2125. 
T.I. acknowledges support from VILLUM Foundation (grant no. VIL37766) and the DNRF Chair program (grant no. DNRF162) by the Danish National
Research Foundation. 
This work is supported in part by JSPS KAKENHI Grant Nos. 20H05850 (C.Y.), 20H05853 (T.H., C.Y.), 21K20367 (Y.K.), 23KK0048 (Y.K.), 24K07027 (T.H., C.Y.), 24KJ1223 (D.S.), and 25K07281 (C.Y.).


# References


