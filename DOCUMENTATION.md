# Ocean Engineering Toolbox (OET) Documentation

The Ocean Engineering Toolbox (OET) is an open-source Modelica library to simulate a single-body, floating, rigid body in uni-directional waves (mono- and polychromatic). The library is under development and currently in pre-release version 0.2. A typical use case for the OET is the simulation of freely-floating wave energy conversion devices.

This documentation applies to the [latest release](https://github.com/ajay-menon-iitkgp/OET_Sys-MoDEL/releases/latest) only. Previous versions of the OET are not maintained and may be unstable on newer versions of Modelica/OpenModelica.

> [!WARNING]
> The state-space BEMIO import statements in `Processor.m` are under construction. Users must manually enter the SS matrices or use the tutorial `hydroCoeff.mat` file when simulating the RM3 float. - \[28 February 2024\]

Quick links:

- [About the OET](#about-the-ocean-engineering-toolbox)
- [Getting Started](#getting-started)
- [Tutorial](#tutorial)

> [!NOTE]
> To use the OET, a user must possess a MATLAB license and have installed the WEC-Sim Code.

## About the Ocean Engineering Toolbox
Modelica is a symbolic programming language developed by the Modelica Association :tm: that represents Cyber-physical Systems (CPS) in the time-domain. Although it has widespread applications in the automobile and aerospace industries, there is no dedicated Modelica Standard Library (MSL) for the ocean engineering community. A reason for this is the difficulty with representing frequency-dependent variables, since Modelica operates solely in the time-domain. Furthermore, the solution to convolution intergals is a significant challenge since a time history of variables cannot be accessed [\[1\]](#publications).

The Cummins equation is selected for the OET as an appropriate frequency-dependent formulation of floating structures. The radiation and excitation forces are represented by convolution integrals and have been succesfully implemented in v0.1 and v0.2 respectively using a State-Space approach. All results are validated against the WEC-Sim time-domain solver developed in SimScape (read about WEC-Sim [here](https://github.com/WEC-Sim/WEC-Sim)).

This GitHub repository consists of the following components:

- `OET.mo` - OET library source code (Modelica)
- `Processor.m` - Pre- and post-processing script (Matlab)
- `DOCUMENTATION.md` - User guide (Markdown)
- `\tutorial` - Sub-directory containing example files
- `\Previous Versions` - Sub-directory with previous OET versions (not maintained)

The OET Modelica library is broken down into packages based on functionality, each with models (i.e., components) to build a simulation. The OET directory is as follows:

- <details open>
  <summary>Ocean Engineering Toolbox</summary>

  - <details>
    <summary>Wave Profile</summary>

    - <details>
      <summary>Regular Wave</summary>
 
      - LinearWave (Monochromatic linear Airy wave)
      
      </details>
    
    - <details>
      <summary>Irregular Wave</summary>
        
      - PiersonMoskowitzWave (Fully-developed sea state)
      
      - BretschneiderWave (Modified PM spectrum for developing sea state)

      - JONSWAPWave (Developing sea state with limited fetch)

      </details>

    </details>
  
  - <details>
    <summary>Structures</summary>

    - RigidBody (Solves 1DOF motion using the Cummins equation)

    </details>
  
  - <details>
    <summary>Internal</summary>
      
    - <details>
      <summary>Functions</summary>
 
      - WaveNumber (Wave number iterations from frequency and depth)
      
      - randomNumberGen (Random numbers through XOR shift generator)

      - frequencySelector (Select wave frequencies from a range)

      - spectrumGenerator_PM (Generate Pierson Moskowitz spectrum for frequency components)

      - spectrumGenerator_BRT (Generate Brettschneider spectrum for frequency components)
        
      - spectrumGenerator_JONSWAP (Generate JONSWAP spectrum for frequency components)

      </details>
    
    - <details>
      <summary>Connectors</summary>
        
      - WaveOutConn (Output transfer wave elevation and excitation force)
      
      - WaveInConn (Input transfer wave elevation and excitation force)

      - DataCollector (Transfer 'Rigid Body' dynamics and forces)

      </details>
    
    - TestDevelopment (Developer component to test all models, functions, connectors)
    
    </details>
  
  - <details>
    <summary>Tutorial</summary>
      
    - Sample 1 (Example model to simulate a rigid body in regular waves)
      
    - Sample 2 (Example model to simulate a rigid body in irregular waves)
    
    </details>
  
  - <details>
    <summary>Simulations</summary>
      
    - (Directory for users to build custom simulation models)
    
    </details>
    
  </details>

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

#### OET Files

The OET consists of the following core files available under release v0.2:

1. `OET.mo`
2. `Processor.m`

Both files must be downloaded locally to run the OET.

1. `OET.mo` contains the library source code and is a Modelica file. Users must import this library into OpenModelica Editor (OMEdit) to build simulations. The library is split into core files, a directory for tutorials, and a directory for custom simulations.
2. `Processor.m` is a MATLAB file that runs pre-processing and post-processing scripts for the OET [input files](#input-file). To run BEMIO, the user must already have installed and setup WEC-Sim. At present, direct processing of output from BEM tools (WAMIT, Nemoh, Capytaine) are not supported. 
    - The `preProcessor` section reads an HDF file generated by WEC-Sim's BEMIO tool. This HDF file contains the hydrostatic and hydrodynamic coefficients for the simulation geometry. The `preProcessor` then exports the required data into a MATLAB structure that can be read by Modelica.
    - The `postProcessor` section reads the output files generated by the OET, runs WEC-Sim, and plots a comparison of the excitation force and heave dynamic response between both tools. If the user does not want to compare results against WEC-Sim, do not run the `postProcessor` section.

> [!NOTE]
> In place of `OET.mo`, you may find an equivalent `OET-vX-Y.mo` Modelica file in the repository, where X.Y is the release number. For example, `OET-v0-2.mo` is the code for version 0.2. For simplicity, this documentation henceforth refers to the library source code as `OET.mo` irrespective of the version.

> [!CAUTION]
> Always execute `Processor.m` using MATLAB's `Run Section` or `Run and Advance` option, NOT `Run`. This is because the `preProcessor` and `postProcessor` scripts are different sections in the same file.

#### Input File

To simulate a geometry, the user must supply the BEMIO output file in the Hierarchical Data Format (HDF) type. Since Modelica is not directly able to access HDF files, the `preProcessor` script is required to convert it into a MATLAB structure (saved as `hydroCoeff.mat`) that can be imported by the OET. Users are directed to the BEMIO [documentation]() for a more detailed tutorial on generated the HDF file.

The `preProcessor` script extracts and writes the following parameters into the `hydroCoeff.mat` file (quantities in *SI* units):

Parameter | Variable Name |Unit
--- | :---: | :---:
Total body mass | `m33` | $kg$
Added mass (infinite frequency) | `Ainf33` | $kg$
Hydrostatic stiffness | `Khs33` | $N/m$
State Matrix | `A` | $N.s/m$
Input Matrix | `B` | $N.s/m$
Output Matrix | `C` | $N.s/m$
Feed-through Matrix | `D` | $N.s/m$
Wave Excitation Coefficient Frequencies | `w` | $rad/s$
Wave excitation coefficient - real component | `FexcRe` |
Wave excitation coefficient - imaginary component | `FexcIm` |

The OET does not currently support visualization capabilities. Hence, a geometry file is not required.

> [!NOTE]
> BEMIO can process the output files from multiple diffraction analysis tools including the open-source [Nemoh](https://lheea.gitlab.io/Nemoh/index.html) and [Capytaine](https://capytaine.github.io/stable/) software.

#### Step-by-Step Guide

1. Download `OET.mo` and `Processor.m` locally.
2. For a chosen geometry, run a diffraction analysis/BEM tool that is compatible with BEMIO (Capytaine, Nemoh, Aqwa, or WAMIT).
3. Run BEMIO on the BEM output of a chosen geometry to generate an HDF file of the hydrostatic and hydrodynamics coefficients/parameters.
4. Run only the `preProcessor` section of the `Processor.m` file. The required parameters are saved to a structure `hydroCoeff.mat` in the user's current MATLAB directory.
5. Startup OpenModelica and import the OET library source code (`OET.mo`).
6. In the `OceanEngineeringToolbox/Simulations` directory, the user can perform a custom simulation by creating a model. Building a simulation consists of connecting a chosen wave excitation component (regular or from a choice of three irregular spectra) to the `RigidBody` model in the `\Structures` OET sub-directory. To do so, create a new model and populate it with instances of the required components. Then, connect the data ports (or *connectors*) of both components using the `connect()` function. The user must input the full path for the `hydroCoeff.mat` file in the parameter `String filePath`.
7. The specifications of a component can be changed within the model by explicitly defining any variables marked as a 'parameter' in the component source code. For instance, the significant height (`Hs`), number of wave frequency components (`n_omega`), and the ramp time (`Trmp`) are three parameters in the irregular wave profile package (the regular wave package sets `n_omega = 1`). A default value for each parameter is provided in the component source code, but a user-defined value can be specified when creating an instance of the wave profile component. The following examples illustrate this and the previous step, and users can append other parameters in a similar fashion:
    ```modelica
    package Tutorial

        /* Model to simulate a body in regular waves */
        model sample1
            parameter String filePath = "F:/.../hydroCoeff.mat";          /* Edit this parameter to full path of hydroCoeff.mat */
            OceanEngineeringToolbox.WaveProfile.RegularWave.LinearModel Reg1(fileName = filePath, Hs = 5, Trmp = 50);    /* Edit component parameters if required */
            OceanEngineeringToolbox.Structures.RigidBody Body1(fileName = filePath);        /* Body1 instance of RigidBody component initialized */

            equation
                connect(Reg1.wconn.F_exc, Body1.wconn.F_exc);              /* Components connected */
                annotation(...);                                          /* Edit the annotation to control the simulation parameters */
        end sample1;

        /* Model to simulate a body in irregular waves - Pierson Moskowitz spectrum */
        model sample2
            parameter String filePath = "F:/.../hydroCoeff.mat";
            OceanEngineeringToolbox.WaveProfile.IrregularWave.PiersonMoskowitzWave PM1(fileName = filePath, Hs = 2.5, n_omega = 100, Trmp = 100);    /* Additional 'n_omega' parameter for irregular waves */
            OceanEngineeringToolbox.Structures.RigidBody Body2(fileName = filePath);

            equation
                connect(PM1.wconn.F_exc, Body2.wconn.F_exc);
                annotation(...);
        end sample2;

    end Tutorial;
    ```

> [!IMPORTANT]
> Note, that simulation parameters must stay fixed for the duration of the simulation, but can be modified between simulations. For instance, the value of `Hs` differs between `sample1` and `sample2` (two different simulations). When the user simulates a model, this value must stay constant and cannot be modified mid-simulation.

8. To change the simulation parameters, the user must modify the `annotation(...)` section of each model. The user must use caution when attempting to modify other parameters (such as the solver or integration methods). The following parameters can be edited:

    - `StartTime` - Virtual start time of the simulation (default = 0 seconds)
    - `StopTime` - Virtual simulation duration incremented to `StartTime` (default = 400 seconds)
    - `Interval` - Time step (default = 0.1 seconds)
    - `Tolerance` - Integration tolerance (default 1e-06)

9. Expand the library browser panel (GUI left), right-select the user-defined simulation model (`sample1` in the above example), and select the `simulate` option to build the model.
10. Once the simulation is over, navigate to the 'Plot' tab (GUI bottom-right) and select the variables to view. The following variables are significant to typical users and can be found using the search bar on the variables panel:

Parameter | Variable Name |Unit
--- | :---: | :---:
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
Wave number of wave components | `k` | $m^{-1}$
Wave elevation profile | `SSE` | $m$
Wave excitation coefficient - real component | `ExcCoeffRe` |
Wave excitation coefficient - imaginary component | `ExcCoeffIm` |

11. These output variables can be exported to a CSV file for post-processing or data analysis. For users who wish to compare results with WEC-Sim, the variables must be selected in a specific order and then exported, so that the CSV output preserves this order. The default Modelica export file name `exportedVariables.csv` is used in `Processor.m`. Selecting any variable automatically saves the `time` variable to the first column of the CSV file. Thus, the user must select the variables from `SSE` onwards. The order of these variables (and the order in the CSV file) is as follows:

    - Time (`t`) - Does not need to be selected
    - Wave elevation profile (`SSE`)
    - Wave excitation force (`F_exc`)
    - Wave radiation force (`F_rad`)
    - Heave displacement response (`z`)
    - Heave velocity response (`v_z`)
    - Heave acceleration response (`a_z`)

12. The `postProcessor` script is divided into multiple sections that must be run sequentially. Script 1 generates an `elevationData.mat` file in the current MATLAB directory. This is the WaveClass input for WEC-Sim to simulate using the same wave profile generated by the OET.
13. Script 2 provides important instructions for the user on setting up the WEC-Sim simulation and then executes WEC-Sim in the current directory when run. These instructions are furthered explained in the tutorial below.
14. Script 3 plots the time series comparisons of various aspects of the simulation (forces and responses) between the OET (from the CSV file) and WEC-Sim. The figures generated by this section are:

    - Wave elevation profile (1)
    - Wave heave excitation force (2.1)
    - Heave radiation force (2.2)
    - Heave displacement response (3.1)
    - Heave velocity response (3.2)
    - Heave acceleration response (3.3)
    - Absolute difference - excitation force (4.1)
    - Absolute difference - radiation force (4.2)
    - Absolute difference - heave velocity (4.3)

> [!TIP]
> Developers and OET contributors can test code modifications by running the `Internal\TestDevelopment` model in the OET source rather than having to manually build and simulate multiple models. `TestDevelopment` builds sub-models with each excitation component and the rigid body. Errors and warnings can be identified and addressed by simulating this model.

## Tutorial

This section provides a step-by-step example of using the OET. The Reference Model 3 (**RM3**) device is a WEC developed under the Reference Model Project (RMP) by NREL & Sandia National Labs, and funded by the U.S. DoE. It consists of a float constrained to linearly float along a spar. The relative translation drives a Power Take-Off (PTO) system. This example models the RM3 float, examines its dynamic response in an irregular, uni-directional wave field, and validates the results using WEC-Sim.

All files used in this tutorial are available in the `\tutorial` directory of this repository. The contents are:

- `hydroCoeff.mat` - MATLAB structure with hydro coefficient for RM3 required by the OET
- `exportedVariables.csv` - OET output dataset for the RM3 float simulation
- `wecSimInputFile.m` - WEC-Sim input file to run simulate the RM3 float
- `float.slx` - Simulink model of the RM3 float (to run WEC-Sim)
- `elevationData.mat` - MATLAB structure containing the wave elevation time series generated by `Processor.m`
- `\hydroData` - WEC-Sim sub-directory
    - `rm3.out` - WAMIT output file of RM3 (float and spar)
    - `bemio.m` - BEMIO Code from WEC-Sim
    - `rm3.h5` - BEMIO output HDF5 dataset
- `\geometry` - WEC-Sim sub-directory
    - `float.stl` - Mesh file of the float geometry
- `\output` - WEC-Sim sub-directory
    - `float_matlabWorkspace.mat` - MATLAB structure with WEC-Sim workspace variables (and output)
- `\Validation Images` - Sub-directory with comparison plots between WEC-Sim and the OET
    - `fig(i).fig` - MATLAB figure #i generated by `Processor.m`
    - `fig(i).png` - PNG figure #i generated by `Processor.m`

#### Step 0 - Setup

1. Download the tutorial directory from the repository to the user's local device. The git clone functionality may be used for this purpose.
2. This local directory shall be labelled '$CASE` and must be set as the MATLAB current directory from which the WEC-Sim tool is run.
3. Download `OET.mo` and `Processor.m` to the `$CASE` directory.

#### Step 1 - Diffraction Analysis and Running BEMIO

1. Perform a diffraction analysis of the RM3 float (`float.stl` mesh provided). `rm3.out` is the WAMIT output file for RM3, where 'body1' is the float and 'body2' is the spar.
2. Open MATLAB and run WEC-Sim's BEMIO Code on `rm3.out` within the `$CASE\hydroData` sub-dir. The output is saved as `rm3.h5`.

#### Step 2 - Running the PreProcessor

1. Open the `Processor.m` MATLAB script and verify the curent directory is set to `$CASE`. Note, each section of this script must be sequentially run.
2. In the `Preprocessor Script` section, change the `filename` variable to the full path address of the `rm3.h5` file saved on the user's local device. This variable is used by the `h5read()` function.
3. Run the `PreProcessor Script` section of the file and advance. This generates the `hydroCoeff.mat` data structure containing all OET simulation parameters in the current directory.
4. The user can compare their preProcessing results with the sample `hydroCoeff.mat` file provided in `$CASE` for the RM3 float geometry.

#### Step 3 - Building and Running the OET Simulation

1. In OpenModelica, import the `OET.mo` library file and open it in the 'text edit' mode.
2. This tutorial simulates the RM3 float in irregular waves defined using the single-peak Pierson Moskowitz (PM) spectrum. The tutorial can be simulated by building model `sample1` in the OET's `OceanEngineeringToolbox\SampleSimulation`.
3. In `sample1`, edit the filePath parameter variable to the full path address of the `hydroCoeff.mat` file saved locally.
4. In `sample1`, edit the required wave component parameters (`Hs`, `Trmp`, `n_omega`) for the PM1 instance.
5. In `sample1`, navigate to the annotation section and edit the simulation parameters as required (`StopTime`, `Interval`, etc.).
6. Save these changes to the model. On the library browser panel of the OpenModelica GUI, expand the OET library, open the `SampleSimulation` package, right-click `sample1`, and select the `simulate` option. This executes the tutorial simulation.
7. Once the simulation finishes, open the 'Plot' tab of the GUI to view the results.
8. Navigate to the 'variable browser' pane. Using the variable nomenclature defined in the '[Getting Started](#getting-started)' section, select the output variables specified in the step-by-step guide in the specified order (`SSE`, `F_exc`, `F_rad`, `z`, `v_z`, `a_z`). Note, that `time` is automatically saved to the first column.
9. Select the `Export to CSV` option to save the results. If the user wishes to perform a comparison study in WEC-Sim, this CSV file must be saved in the same directory as the WEC-Sim simulation. The default file name must be used (`exportedVariables.csv`).
10. If the user wishes to revisit the results of the OET, they can save the simulation workspace in OpenModelica.

#### Step 4 - Running the PostProcessor

1. In the `Processor.m` file, navigate to the section titled 'Next Steps'. Follow the instructions to run the same simulation in WEC-Sim and compare the results to the OET.
2. In the `PostProcessor Script - 1` section, ensure that the CSV file specified in the filename argument of `readmatrix()` is present in the working WEC-Sim directory. Verify that `$CASE` is the MATLAB current directory.
3. Run this section to generate and save the wave elevation time series data structure labelled `elevationData.mat` in the `$CASE` directory.
4. Prior to running the `PostProcessor Script - 2` section, the user must setup the WEC-Sim simulation. This entails the following steps:

    - In `wecSimInputFile.m`, uncomment the `elevationImport` wave class and comment out all other wave classes. Ensure that `waves.elevationFile = 'elevationData.mat';` is defined.
    - Verify that the `elevationData.mat` and `exportedVariables.csv` files are in the same `$CASE` directory as the WEC-Sim simulation.
    - The remaining setup involves the standard WEC-Sim workflow - create a Simulink model `float.slx`, create the `\hydroData` sub-dir with the `rm3.h5` HDF5 file, and create the `\geometry` sub-dir with the 'float.stl' mesh file saved in it. The user can refer to the repository's `\tutorial` directory for the example setup.

5. Run the WEC-Sim simulation from the `$CASE` directory by either executing the `wecSim` function on the MATLAB command line or running the `PostProcessor Script - 2` section of `Processor.m`.
6. Once the WEC-Sim simulation finishes, the user can plot a series of comparisons by running the `PostProcessor Script - 3` section. Sample plots are included in the `\tutorial\Validation Images` directory of this repository. Note that the final results may vary since a random wave elevation profile is generated each time. To generate repeatable waves, the user can fix the local and global seeds defined for each wave component in the OET source code.

Feedback, clarifications, and issues may be brought to the attention of the development team using the 'Issues' board of this repository.
