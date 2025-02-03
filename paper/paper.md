---
title: 'COSMOS: A c++ package for numerical relativity specialized for PBH formation'
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
   index: 5
date: 3 February 2025
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Primordial black holes (PBHs) are those generated in the early universe without experience of the form of a star. 
It has been pointed out that PBHs may be candidates for black holes and compact objects of various masses in the universe, 
or may be a major component of dark matter, and are attracting attention.
In the standard formation process, PBHs are formed from super-horizon primordial fluctuation with 
a non-linearly large amplitude. 
In order to follow the whole non-linear gravitational dynamics, one has to rely on numerical relativity 
solving Einstein equations numerically. 
Furthermore, since there is a hierarchy between the size of collapsing region and 
cosmological expansion scale, an efficient resolution refinement procedure is needed.  
`COSMOS` [@Yoo:2013yea; @Okawa:2014nda] and `COSMOS-S` [@Yoo:2021fxs] provide simple tools for simulation of PBH formation. 
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

`COSMOS` and `COSMOS-S` are C++ packages for solving Einstein equations in 3+1 dimension and spherical symmetry, respectively.  
It is originally translated from SACRA code [@Yamamoto:2008js] into C++. 
Perfect fluid with linear equation of states and massless scalar field are implemented as matter fields. 
In order to resolve the collapsing region, non-Cartesian scale-up coordinates [@Yoo:2018pda] and 
a fixed mesh-refinement procedure are implemented. 
OpenMP package is used for the parallelization. 
No other packages are not required, but users are supposed to understand the source code to some extent and use it by modifying it themselves. 


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

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References