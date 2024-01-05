## Ocean Engineering Toolbox (OET) - Modelica Library
The Ocean Engineering Toolbox (OET) is a Modelica library to simulate a single-body, floating Wave Energy Converter (WEC) in uni-directional, mono- and polychromatic waves.

This is a work in progress. Current pre-release version: OET v0.2 (2024).

#### Overview

Modelica is a symbolic programming language used to represent Cyber-physical Systems (CPS) in the time-domain. Although it has widespread applications in the automobile and aerospace industries, there is no dedicated Modelica Standard Library (MSL) for the ocean engineering community. A significant reason is the challenge with representing frequency-dependent variables, since Modelica operates solely in the time-domain. Furthermore, the solution to convolution intergals is a significant challenge since a time history of variables cannot be accessed.

Literature demonstrates the critical nature of fequency-dependent models, and the Cummins equation is selected for the OET as an appropriate formulation. The radiation and excitation forces are represented by convolution integrals and have been succesfully implemented in v0.1 and v0.2 respectively. All results are validated against the WEC-Sim time-domain solver developed in SimScape (read about WEC-Sim [here](https://github.com/WEC-Sim/WEC-Sim)).

#### Development Team
The OET is a product of research at the Sys-MoDEL research group, University of New Brunswick - Canada. The researchers are:

- Ajay Menon
- Ali Shahbaz Haider
- Kush Bubbar

#### Citing This Work
Publications relevant to this tool are:

- A. Menon, A. S. Haider, and K. Bubbar. "On Implementing Cummins Equation to Represent Accurate Wave Radiation Forces in Modelica." in *Proceedings of the 42nd International Conference on Offshore Mechanics and Arctic Engineering*, vol. 86908, p. V008T09A068. American Society of Mechanical Engineers, 2023.

#### Acknowledgements
The team would like to thank MITACS Canada for funding the collaborative research with the Sys-MoDEL research team at the University of New Brunswick via the Globalink 
Research Internship program in 2022 and 2024. The developers also thank Maplesoft™ for their technical applications support to this work.
