# Ocean Engineering Toolbox (OET) Documentation

The Ocean Engineering Toolbox (OET) is a Modelica library to simulate a single-body, floating, rigid body in uni-directional waves (mono- and polychromatic). The library is under development and currently in pre-release version 0.2. A typical use case for the OET is the simulation of freely-floating wave energy conversion devices. This documentation covers the installation and usage of v0.2 when running the OpenModelica (OMEdit) solver engine.

Quick links:

- [About the OET](#about-the-ocean-engineering-toolbox)
- [Getting Started](#getting-started)
- [Tutorial](#tutorial)

> [!NOTE]
> The documentation applies to the [latest release](https://github.com/ajay-menon-iitkgp/OET_Sys-MoDEL/releases/latest) only. Previous versions of the OET are not maintained and may be unstable on newer versions of Modelica/OpenModelica.

## About the Ocean Engineering Toolbox
Modelica is a symbolic programming language developed by the Modelica Association that represents Cyber-physical Systems (CPS) in the time-domain. Although it has widespread applications in the automobile and aerospace industries, there is no dedicated Modelica Standard Library (MSL) for the ocean engineering community. A reason for this is the difficulty with representing frequency-dependent variables, since Modelica operates solely in the time-domain. Furthermore, the solution to convolution intergals is a significant challenge since a time history of variables cannot be accessed [\[1\]](#publications).

The Cummins equation is selected for the OET as an appropriate frequency-dependent formulation of floating structures. The radiation and excitation forces are represented by convolution integrals and have been succesfully implemented in v0.1 and v0.2 respectively using a State-Space approach. All results are validated against the WEC-Sim time-domain solver developed in SimScape (read about WEC-Sim [here](https://github.com/WEC-Sim/WEC-Sim)).

#### Development Team
The OET is a research product from the Sys-MoDEL group at the University of New Brunswick (Fredericton), Canada. Current members of the development team are:

- Ajay Menon
- Ali Shahbaz Haider
- Kush Bubbar

#### Publications
Publications on the development of the OET are:

\[1\] A. Menon, A. S. Haider, and K. Bubbar. "On Implementing Cummins Equation to Represent Accurate Wave Radiation Forces in Modelica." in *Proceedings of the 42nd International Conference on Offshore Mechanics and Arctic Engineering*, OMAE 2023, Melbourne, Australia, 2023.

\[2\] A. Menon, A. S. Haider, and K. Bubbar. "Advancing the Modelica Ocean Engineering Toolbox with the Capability to Generate Accurate Wave Excitation Forces." [ACCEPTED] *43rd International Conference on Offshore Mechanics and Arctic Engineering*, OMAE 2024, Singapore, 2024.

## Getting Started

This section outlines the workflow for the OET and a step-by-step guide on running simulations. Note, that v0.2 of the OET only simulates single body structures that are freely-floating and subject to uni-directional waves (regular or irregular).

#### Setting Up OpenModelica

As a scientific programming language, Modelica requires a Modelica Solver Engine (MSE) to build and execute simulations. OpenModelica is an **open-source** MSE and the OET is developed to run primarily on it. By modifying the code annotations, other alternative engines such as MapleSim, Modelon, and Dymola can also be used. OpenModelica is easy to install and users are directed to download the latest version from its [website](https://openmodelica.org). The Modelica [Documentation](https://openmodelica.org/doc/OpenModelicaUsersGuide/1.22/) guide provides a detailed list of libraries, functionalities, and tutorials for new users to learn the language.

To build a model from library components, users must define an instance of each component they will use. This can be done by referencing the relative address of the component in the library browser window. Examples of the instance definition used in the OET are:

```
Instance definition examples here
```

Models require the components to be assembled, and this is achieved using the `connect()` function, as follows:

```
Examples of connect() function here
```

> [!IMPORTANT]
> Users require a MATLAB license to run the `Processor.m` script and may require additional toolboxes if they wish to run the WEC-Sim toolbox in SimScape for comparison/validation purposes.

#### OET Files
The OET consists of the following core files available under the relevant release:

1. `OET.mo` (or an equivalent OET-vX-Y.mo file, where X.Y is the latest release. For example, `OET-v0-2.mo` is the code for version 0.2. For simplicity, the documentation references this file as `OET.mo` henceforth.)
2. `Processor.m`

Both files must be downloaded locally to run the OET.

1. `OET.mo` contains the library source code and is a Modelica file. Users must import this library into OpenModelica Editor (OMEdit) to build simulations. The library is split into core files, a directory for tutorials, and a directory for custom simulations.
2. `Processor.m` is a MATLAB file that runs pre-processing and post-processing scripts for the OET [input files](#required-input-file).
    - The `preProcessor` section reads an HDF file generated by WEC-Sim's BEMIO tool. This HDF file contains the hydrostatic and hydrodynamic coefficients for the simulation geometry. The `preProcessor` then exports the required data into a MATLAB structure that can be read by Modelica.
    - The `postProcessor` section reads the output files generated by the OET, runs WEC-Sim, and plots a comparison of the excitation force and heave dynamic response between both tools. To run `postProcessor.m`, the user must already have installed and setup WEC-Sim. At present, other tools are not supported. If the user does not want to compare results against WEC-Sim, do not run the `postProcessor` section.

> [!CAUTION]
> Always execute `Processor.m` using MATLAB's `Run Section and Advance` option, not `Run`. This is because the `preProcessor` and `postProcessor` scripts are different sections in the same file.

#### Input File

To simulate a geometry, the user must supply the BEMIO output file in the Hierarchical Data Format (HDF) type. Since Modelica is not directly able to access HDF files, the `preProcessor` script is required to convert it into a MATLAB structure (saved as `hydroCoeff.mat`) that can be imported by the OET. Users are directed to the BEMIO [documentation]() for a more detailed tutorial on generated the HDF file.

The `preProcessor` script extracts and writes the following parameters into the `hydroCoeff.mat` file (quantities in *SI* units):

Parameter | Variable Name
--- | ---
Total body mass | `M`
Added mass (infinite frequency) | `Ainf`
Hydrostatic stiffness | `C`
State Matrix \[**A**\] | `A1`
Input Matrix \[**A**\] | `B1`
Output Matrix \[**C**\] | `C1`
Feedthrough/feedforward Matrix \[**D**\] | `D1`
Wave Excitation Coefficient Frequencies \[rad/s\] | `w`
Wave excitation coefficient - real component | `FexcRe`
Wave excitation coefficient - imaginary component | `FexcIm`

The OET does not currently support visualization capabilities. Hence, a geometry file is not required.

> [!NOTE]
> BEMIO can process the output files from multiple diffraction analysis tools including the open-source [Nemoh]() and [Capytaine]() software.

#### Step-by-Step Guide

1. Download `OET.mo` and `Processor.m`.
2. Run BEMIO on the diffraction analysis output of a chosen geometry to generate an HDF file of the hydrostatic and hydrodynamics coefficients/parameters.
3. Run only the `preProcessor` section of the `Processor.m` file. The required parameters are saved to a structure `hydroCoeff.mat`.
4. Startup OpenModelica and import the OET as a library file.
5. Find all instances of `Modelica.Utilities.Streams.readRealMatrix("file-address/hydroCoeff.mat", "hydroCoeff.variable", 1, size1, size2)` by searching for the function `readRealMatrix`. Replace `file-address` with the the full path for `hydroCoeff.mat`. For variable names, refer to [this table](#required-input-file) and follow the same nomenclature. Refer to the HDF file to identify numeric values for `size1` and `size` of the variable. For a scalar variable (not an array/list), set both sizes to `1`.

> [!TIP]
> For example, `Modelica.Utilities.Streams.readRealMatrix("C:/Users/OET/Example/hydroCoeff.mat", hydroCoeff.FexcRe, 1, 260)`.

7. Create a new model within the code by following the example model under package `SampleSimulations`. The user must connect a wave excitation package (regular, PM, Bretschneider, Jonswap) with the WEC model.
8. > [!NOTE]
   > Provide more information on writing code for a new model and give some examples.
9. 
10. Run the simulation by right-clicking the user-defined model on the libraries browser panel (located to the left) and selecting 'Build'.
11. Once the simulation has been executed, head to the 'Plot' tab and select the variables to view.
12. > [!NOTE]
    > Provide more information on exporting, the variables to export, nomenclature of variables, etc.
13. To compare results against WEC-Sim, the `postProcessor` section of `Processor.m` must be run. Prior to this, export the ``, ``, ``, and `` variables to a CSV file from the OET variables browser.
14. Run the `postProcessor` sub-sections iteratively to first execute WEC-Sim and then plot the comparisons.

> [!IMPORTANT]
> To compare the results with WEC-Sim, the time series of the wave elevation generated by the OET must be imported into WEC-Sim. Instructions to do this and setup the WEC-Sim input files are provided in `Processor.m`.

## Tutorial
This section provides users with an exmaple to refer to when building and executing simulations in the OET.

#### Geometry and Diffraction Analysis

#### Running the PreProcessor

#### Building and Running the OET Simulation

#### Running the PostProcessor
