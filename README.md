## Ocean Engineering Toolbox (OET) - Modelica Library
The Ocean Engineering Toolbox (OET) is a Modelica library to simulate a single-body, floating Wave Energy Converter (WEC) in uni-directional, mono- and polychromatic waves.

This is a work in progress. Current pre-release version: OET v0.2 (2024).
[Documentation](https://ajay-menon-iitkgp.github.io/OET_Sys-MoDEL/) under construction. User guide available [here](https://github.com/ajay-menon-iitkgp/OET_Sys-MoDEL/wiki/Ocean-Engineering-Toolbox-(OET)-%E2%80%90-Modelica-Hydrodynamics-Library).

#### Overview

Modelica is a symbolic programming language used to represent Cyber-physical Systems (CPS) in the time-domain. Although it has widespread applications in the automobile and aerospace industries, there is no dedicated Modelica Standard Library (MSL) for the ocean engineering community. A reason for this is the difficulty with representing frequency-dependent variables, since Modelica operates solely in the time-domain. Furthermore, the solution to convolution intergals is a significant challenge since a time history of variables cannot be accessed.

Literature demonstrates the requirement for fequency-dependent models, and the Cummins equation is selected for the OET as an appropriate formulation. The radiation and excitation forces are represented by convolution integrals and have been succesfully implemented in v0.1 and v0.2 respectively using a State-Space approach. All results are validated against the WEC-Sim time-domain solver developed in SimScape (read about WEC-Sim [here](https://github.com/WEC-Sim/WEC-Sim)).

#### Development Team
The OET is a research product from the Sys-MoDEL group at the University of New Brunswick (Fredericton), Canada. Current members of the development team are:

- Ajay Menon
- Ali Shahbaz Haider
- Kush Bubbar

#### Citing This Work
Publications on the development of the OET are:

- A. Menon, A. S. Haider, and K. Bubbar. "On Implementing Cummins Equation to Represent Accurate Wave Radiation Forces in Modelica." in *Proceedings of the 42nd International Conference on Offshore Mechanics and Arctic Engineering*, OMAE 2023, Melbourne, Australia, 2023.
- A. Menon, A. S. Haider, and K. Bubbar. "Advancing the Modelica Ocean Engineering Toolbox with the Capability to Generate Accurate Wave Excitation Forces." [ACCEPTED] *43rd International Conference on Offshore Mechanics and Arctic Engineering*, OMAE 2024, Singapore, 2024.
