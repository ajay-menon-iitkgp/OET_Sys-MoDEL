# Ocean Engineering Toolbox (OET) Documentation

The Ocean Engineering Toolbox (OET) is an open-source Modelica library to simulate a single-body, floating, rigid body in uni-directional waves (mono- and polychromatic). The library is under development and currently in pre-release version 0.2. A typical use case for the OET is the simulation of freely-floating wave energy conversion devices. This documentation covers the installation and usage of v0.2 when running the OpenModelica (OMEdit) solver engine.

Quick links:

- [About the OET](#about-the-ocean-engineering-toolbox)
- [Getting Started](#getting-started)
- [Tutorial](#tutorial)
- [Advanced Concepts](#advanced-concepts)

> [!NOTE]
> The documentation applies to the [latest release](https://github.com/ajay-menon-iitkgp/OET_Sys-MoDEL/releases/latest) only. Previous versions of the OET are not maintained and may be unstable on newer versions of Modelica/OpenModelica.

## About the Ocean Engineering Toolbox
Modelica is a symbolic programming language developed by the Modelica Association :tm: that represents Cyber-physical Systems (CPS) in the time-domain. Although it has widespread applications in the automobile and aerospace industries, there is no dedicated Modelica Standard Library (MSL) for the ocean engineering community. A reason for this is the difficulty with representing frequency-dependent variables, since Modelica operates solely in the time-domain. Furthermore, the solution to convolution intergals is a significant challenge since a time history of variables cannot be accessed [\[1\]](#publications).

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

#### OET Files

The OET consists of the following core files available under release v0.2:

1. `OET.mo`
2. `Processor.m`

Both files must be downloaded locally to run the OET.

1. `OET.mo` contains the library source code and is a Modelica file. Users must import this library into OpenModelica Editor (OMEdit) to build simulations. The library is split into core files, a directory for tutorials, and a directory for custom simulations.
2. `Processor.m` is a MATLAB file that runs pre-processing and post-processing scripts for the OET [input files](#input-file).
    - The `preProcessor` section reads an HDF file generated by WEC-Sim's BEMIO tool. This HDF file contains the hydrostatic and hydrodynamic coefficients for the simulation geometry. The `preProcessor` then exports the required data into a MATLAB structure that can be read by Modelica.
    - The `postProcessor` section reads the output files generated by the OET, runs WEC-Sim, and plots a comparison of the excitation force and heave dynamic response between both tools. To run `postProcessor.m`, the user must already have installed and setup WEC-Sim. At present, other tools are not supported. If the user does not want to compare results against WEC-Sim, do not run the `postProcessor` section.

> [!NOTE]
> In place of `OET.mo`, you may find an equivalent `OET-vX-Y.mo` Modelica file, where X.Y is the release number. For example, `OET-v0-2.mo` is the code for version 0.2. For simplicity, this documentation henceforth refers to the library source code as `OET.mo` irrespective of the version.

> [!CAUTION]
> Always execute `Processor.m` using MATLAB's `Run Section and Advance` option, not `Run`. This is because the `preProcessor` and `postProcessor` scripts are different sections in the same file.

> [!IMPORTANT]
> Users require a MATLAB license to run the `Processor.m` script and may require additional toolboxes if they wish to run the WEC-Sim toolbox in SimScape for comparison/validation purposes.

#### Input File

To simulate a geometry, the user must supply the BEMIO output file in the Hierarchical Data Format (HDF) type. Since Modelica is not directly able to access HDF files, the `preProcessor` script is required to convert it into a MATLAB structure (saved as `hydroCoeff.mat`) that can be imported by the OET. Users are directed to the BEMIO [documentation]() for a more detailed tutorial on generated the HDF file.

The `preProcessor` script extracts and writes the following parameters into the `hydroCoeff.mat` file (quantities in *SI* units):

Parameter | Variable Name |Unit
--- | :---: | :---:
Total body mass | `M` | $kg$
Added mass (infinite frequency) | `Ainf` | $kg$
Hydrostatic stiffness | `C` | $N/m$
State Matrix | `A1` | $N.s/m$
Input Matrix | `B1` | $N.s/m$
Output Matrix | `C1` | $N.s/m$
Feed-through / feed-forward Matrix | `D1` | $N.s/m$
Wave Excitation Coefficient Frequencies | `w` | $rad/s$
Wave excitation coefficient - real component | `FexcRe` |
Wave excitation coefficient - imaginary component | `FexcIm` |

The OET does not currently support visualization capabilities. Hence, a geometry file is not required.

> [!NOTE]
> BEMIO can process the output files from multiple diffraction analysis tools including the open-source [Nemoh](https://lheea.gitlab.io/Nemoh/index.html) and [Capytaine](https://capytaine.github.io/stable/) software.

#### Step-by-Step Guide

1. Download `OET.mo` and `Processor.m`.
2. For a chosen geometry, run a diffraction analysis/BEM tool that is compatible with BEMIO (Capytaine, Nemoh, etc.).
3. Run BEMIO on the BEM output of a chosen geometry to generate an HDF file of the hydrostatic and hydrodynamics coefficients/parameters.
4. Run only the `preProcessor` section of the `Processor.m` file. The required parameters are saved to a structure `hydroCoeff.mat` in the user's MATLAB current directory.
5. Startup OpenModelica and import the OET as a library file.
6. Find all instances of `Modelica.Utilities.Streams.readRealMatrix("file-address/hydroCoeff.mat", "hydroCoeff.variable", size1, size2)` by searching for the function `readRealMatrix` in the library source code. Replace `file-address` with the the full path for `hydroCoeff.mat`. For variable names, refer to [this table](#input-file) and follow the same nomenclature. Refer to the HDF file to identify numeric values for `size1` and `size` of the variable. For a scalar variable (not an array/list), set both sizes to `1`. For example,
   ```modelica
   Real Fexc_Re[1, :] = Modelica.Utilities.Streams.readRealMatrix("C:/hydroCoeff.mat", hydroCoeff.FexcRe, 1, 260)
   ```

7. Under the `OceanEngineeringToolbox/Simulations` directory, the user can create a model by connecting a chosen wave excitation component (regular or from a choice of three irregular spectra) to the WEC component. To do so, create a new model and populate it with instances of the required components. Then, connect the data ports (or *connectors*) of both components using the `connect()` function. The following examples illustrate this process to create simulation models,
    ```modelica
    package Simulations

        /* Model to simulate a body in regular waves */
        model sample1
            OceanEngineeringToolbox.WaveProfile.RegularWave.AiryWave Reg1;
            OceanEngineeringToolbox.WEC WEC1;

            equation
              connect(Reg1.wconn.F_exc, WEC1.wconn.F_exc);
        end sample1;

        /* Model to simulate a body in irregular waves - Pierson Moskowitz spectrum */
        model sample2
            OceanEngineeringToolbox.WaveProfile.IrregularWave.PiersonMoskowitzWave PM1;
            OceanEngineeringToolbox.WEC WEC1;

            equation
              connect(PM1.wconn.F_exc, WEC1.wconn.F_exc);
        end sample2;

    end Simulations;
    ```

8.  With the components assembled, expand the library browser panel (GUI left) and right-select the simulation model. Click build and wait for the simulation to complete.
9.  Once the simulation is over, navigate to the 'Plot' tab (GUI bottom-right) and select the variables to view. Owing to the conversion of certain discrete variables to a continuous formulation, redundant intermediate variables may have been created. The following variables are of concern to users:

Parameter | Variable Name |Unit
--- | :---: | :---:
Simulation time | `t` | $s$
Wave ramp time | `Trmp` | $s$
Wave excitation force | `F_exc` | $N$
Wave radiation force | `F_rad` | $N$
Heave displacement response | `z` | $m$
Heave velocity response | `v_z` | $m/s$
Heave acceleration response | `a_z` | $m/s^2$
Number of wave components | `n_omega` | $-$
Wave frequency components | `omega` | $rad/s$
Wave energy spectrum value | `S` | $m^2$
Amplitude of wave components | `zeta` | $m$
Ramped amplitude of wave components | `zeta_rmp` | $m$
Wave number of wave components | `k` | $m^{-1}$
Wave elevation profile, ramped | `SSE` | $m$
Wave excitation coefficient - real component | `ExcCoeffRe` |
Wave excitation coefficient - imaginary component | `ExcCoeffIm` |

10. These output variables can be exported to a CSV file for post-processing or data analysis. For users who wish to compare results with WEC-Sim, the variables must be selected in a specific order and then exported, so that the CSV output preserves this order. The default Modelica export file name `exportedVariables.csv` is used in `Processor.m`. The order of these variables are as follows:

    - Time (`t`)
    - Wave elevation profile, ramped (`SSE`)
    - Wave excitation force (`F_exc`)
    - Wave radiation force (`F_rad`)
    - Heave displacement response (`z`)
    - Heave velocity response (`v_z`)
    - Heave acceleration response (`a_z`)
    - Wave frequency components (`omega`)
    - Energy spectrum value (`S`)
    - Wave excitation coefficient - real component (`ExcCoeffRe`)
    - Wave excitation coefficient - imaginary component (`ExcCoeffIm`)

11. The `postProcessor` script is divided into multiple sections that must be run sequentially. Script 1 generates an `elevationData.mat` file in the current MATLAB directory. This is the WaveClass input for WEC-Sim to simulate using the same wave profile generated by the OET.
12. Script 2 provides important instructions for the user on setting up the WEC-Sim simulation and then executes WEC-Sim in the current directory when run.
13. Script 3 plots the time series comparisons of the heave displacement response and the wave excitation force between the OET (from the CSV file) and WEC-Sim.

## Tutorial
This section provides users with an exmaple simulation in the OET.

#### Geometry and Diffraction Analysis

#### Running the PreProcessor

#### Building and Running the OET Simulation

#### Running the PostProcessor

## Advanced Concepts

Coming soon...