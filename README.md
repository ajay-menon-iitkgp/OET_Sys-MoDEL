# Ocean Engineering Toolbox (OET) - Modelica Library
The Ocean Engineering Toolbox (OET) in Modelica simulates a single-body, floating Wave Energy Converter (WEC).

## Introduction and Background
Modelica is a symbolic programming language used to represent a Cyber-physical System (CPS). Although it has widespread applications in the automobile and aerospace industries, use by the ocean engineering community is negligible. One possible reason is the lack of a Modelica Standard Library that accurately simulates the dynamics of offshore structures.

Wave Energy Converters (WECs) are energy-absorption devices that convert wave energy to an equivalent electrical form, and are an important element in the global renewable energy industry. WECs comprise off the following subsystems:
- Hydrodynamic Subsystem:     Wave-absorption mechanism (captor body).
- Power Take-Off Subsystem:   Converts the motion of the WEC into an electrical energy.
- Reaction Subsystem:         Transfers forces from the mooring and PTO reactions.
- Control Subsystem:          Intelligently "controls" the WEC to maximize energy capture by the PTO subsystem.

A simulation toolbox must be capture frequency-dependent forces that act on the WEC and the Cummins equation is commonly used to represent a floating/moored WEC. The radiation and excitation forces are represented by convolution integrals that capture the history of both forces over the simulation duration. However, no suitable method exists to symbolically represent convolution integrals.

This work seeks to develop a proof-of-concept Modelica library that simulates a single-body WEC subject to a preset wave excitation force. Results closely match data from the benchmark: WEC-Sim (an open-source Simscape Multibody tool developed by researchers at NREL and Sandia Labs, and funded by the U.S. Department of Energy).

## Simulation Setup
Download and import the "WECHydrodynamics.mo" library package to a Modelica solver engine (OMEdit, MapleSim). The following components are imported into the GUI workspace:
- Import Excitation Force
- Cummins Plant Model
- Linearized Validation Toolkit

## Version History
### v1.0 - Linearized Hydrodynamic Model
### v2.0 - Frequency-Dependent Radiation Force (Convolution Integral)

## Future Objectives
Development of a component to output the excitation force data from a chosen wave elevation profile. Users will have an option to with select a wave spectrum model (Pierson-Moskowitz, JONSWAP, Brettschneider) or import their own wave elevation profile.

Implementing a wave excitation convolution is more challenging than the radiation convolution integral, primarily because the excitation force is acausal i.e., the profile has to be forecast for a region to accurately output the excitation force. Future objectives include investigating the use of a Model Predictive Control (MPC) scheme or the application of a Kalman filter to forecast the elevation data.

## Development Team
The researchers involved in this work are:
- Ajay Menon, IIT Kharagpur, India
- Ali Shahbaz Haider, University of New Brunswick, Canada
- Kush Bubbar, University of New Brunswick, Canada

## How to Cite This Work?
This work is currently under review by the Offshore Mechanics and Arctice Engineering (OMAE) Conference - 2023, organized by the American Society of Mechanical Engineers (ASME).

## Acknowledgements
The team would like to thank MITACS Canada for funding the collaborative research with the Sys-MoDEL research team at the University of New Brunswick via the Globalink 
Research Internship program. All the authors graciously thank Maplesoft™ for their technical applications support to this work.
