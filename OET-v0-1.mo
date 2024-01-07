/*  Modelica Ocean Engineering Toolbox v0.1
    Copyright (C) 2022  Ajay Menon

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    Your copy of the GNU General Public License can be viewed
    here <https://www.gnu.org/licenses/>.
*/

/* Library structure:
   OET (LIBRARY)
   |->  Excitation Force (PACKAGE)
   |    |-> Import Excitation Force (MODEL) [!- Import user-defined wave excitation force]
   |
   |->  Plant Model (MODEL)
   |    |-> TF2DE (MODEL) [!- If the user has the transfer function ODE representation]
   |    |-> SS2DE (MODEL) [!- Default method of solving Cummins' equation]
   |
   |->  Reference Models (PACKAGE)
   |    |-> Ref Data 1 (MODEL) [!- Imported from Matlab struct]
   |    |-> Frequency Independent Comparison (MODEL) [!- Only for comparison purposes]
   |
   |->  Functions (PACKAGE) - DEPRECATED
   |
   |->  Connectors (PACKAGE)
        |-> Data Collector (CONNECTOR) [!- Transfer non-flange data]
*/

/* Important notes:
   
   1. To change the library name, replace all instances of "OET" with your chosen name.
   2. The statement "extends Modelica.Blocks.Icons.Block;" is used to define visibility of the model as a block in the Modelica solver. Removing this statement causes the relevant model to disappear from the library.
   3. This library can be run in any Modelica solver engine, although it is intended for use in MapleSim. Annotations have been prepared accordingly.

    Present capabilities:
    1. The present library is capable of importing data from a struct in MATLAB to define the hydrodynamic coefficients, excitation force, and state-space quadruple
    2. Capable of estimating the heave dynamics and radiation force acting on a single body defined using the state space quadruple
*/

/*  How to use the library:
    1.  Download this source code file to any location locally; it does not need to be in Modelica's main directory.
    2.  Run a BEM Code (Nemoh, WAMIT, Capytaine, AQWA) on the desired single-body geometry and import the results into the HDF format using WEC-Sim's BEMIO tool.
    3.  The following parameters are required from the HDF file to run this library on any geometry:
          (i)   Mass (m33)
          (ii)  Infinite added mass (Ainf33)
          (iii) Hydrostatic restoring coefficient (Khs33)
          (iv)  Wave excitation force (F_wav)
          (v)   Validation dataset: actual velocity (vRef)
          (vi)  Validation dataset: actual radiation force (FRadRef)
          (vii) State Space matrices for the radiation force convolution integral
    4.  The above parameters must be stored in a single MATLAB struct and saved to any location.
    5.  In the Modelica code, identify the code "Modelica.Utilities.Streams" and replace the address with the location of the BEM results .mat file in the specified format. Note, the complete address must be entered, and not merely the file's name.
    6.  Import the library into MapleSim using "File -> Modelica -> Import library and documentation".
    7.  Connect the 'ImportExcitationForce' model to the left interface of the 'CumminsPlant' model.
    8.  If you wish to validate the solver's results, connect the right Flange interface of 'CumminsPlant' to the 'RefData' model.
    9.  Set the parameter: 'simulation time' in the solver toolbar and on the parameter panels for components 'ExctitationForce' and 'RefData' (variable: 'Tsim').
    10. MAPLESIM ALERT: Attach probes to monitor the required variables and run the simulation. To view radiation force, attach a probe to the component and NOT the flange. From the flange, only mechanical properties can be observed- displacement (s), velocity (v), acceleration (a), force on the flange (f).
*/

/* Source Code */
package OET
/* Library name: OET */

  package ExcitationForce
  /* Package to process the wave excitation force */
  
    model ImportExcitationForce
      extends Modelica.Blocks.Icons.Block;
      /* Model to import a user-defined wave excitation force as a MATLAB struct
         "Utilities.Streams" library imports MATLAB matrix data and requires full file address
         Converts the discrete import samples to a continuous excitation force */
      
      /* Define output interface for the continuous wave excitation force variable */
      Modelica.Blocks.Interfaces.RealOutput F_exc2 "Continuous, real excitation force signal" annotation(
        Placement(transformation(extent = {{100,-10},{120,10}})));

      /* Component parameters */
      parameter Modelica.SIunits.Time Ts = 0.1 "Discrete sampling time interval in seconds";
      parameter Modelica.SIunits.Time Tsim = 30 "Simulation duration in seconds";

      /* Component variables */
      Real F_exc_discrete2[:,1] = Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.F_wav",5001,1) "5001-samples discrete excitation force at 'Ts' sampling duration";
      Real F_exc_discrete[(Tsim/Ts) + 1,1] = F_exc_discrete2[1:(Tsim/Ts)+1,:] "Discrete excitation force subset based on desired simulation duration 'Tsim'";
      
      equation
        /* Convert the discrete excitation force to a continuous representation
           Loop through the vector and hold the value between the sampling time */
        for k in 1:size(F_exc_discrete,1) loop
          when time > Ts*(k-1)-Ts/2 and time <= Ts*k-Ts/2 then
            F_exc2 = F_exc_discrete[k,1];
          end when;
        end for;
      
    end ImportExcitationForce;
  end ExcitationForce;

  package RadForceModels
  /* Package with components representing different approaches to solve the radiation force convolution integral symbolically
     Documented, but users can opt for the CumminsPlant component instead for a better defined component */
  
    model TF2DEModel
      extends Modelica.Blocks.Icons.Block;
      /* Transfer function implementation of Cummins equation to represent floating bodies
         Transfer function representation of the radiation force convolution integral
         Using this component requires the TF's differential equation to be coded-in */
      
      /* Input interface for the wave excitation force variable */
      Modelica.Blocks.Interfaces.RealInput F_wav "Input flange transferring wave excitation force" annotation(
        Placement(transformation(extent = {{-90,-10},{-110,10}})));
  
      /* Output interface for the body heave dynamics and radiation force */
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b "Output flange transferring force and plant motion" annotation(
        Placement(transformation(extent = {{90,-10},{110,10}})));
  
      /* Component parameters */
      parameter Modelica.SIunits.Mass M = 727010 "Total mass of the body (including ballast)";
      parameter Modelica.SIunits.Mass Ainf = 1232838 "Added mass at maximum (cut-off) frequency";
      parameter Modelica.SIunits.TranslationalSpringConstant C = 2800951.2 "Hydrostatic stiffness";
      /* Replace below parameters for the transfer function coefficients */
      parameter Real a0 = 1 "Coefficient a0 of second-order radiation force derivative";
      parameter Real a1 = 0.9366 "Coefficient a1 of first-order radiation force derivative";
      parameter Real a2 = 1.011 "Coefficient a2 of zeroth-order radiation force derivative";
      parameter Real b0 = 683237 "Coefficient b0 of first-order velocity derivative";
      parameter Real b1 = 54478 "Coefficient b1 of zeroth-order velocity derivative";
  
      /* Component variables */
      Modelica.SIunits.Length z "Heave displacement";
      Modelica.SIunits.Velocity v_z "Heave velocity";
      Modelica.SIunits.Acceleration a_z "Heave acceleration";
      Modelica.SIunits.Force F_rad "Radiation Force";
      
      initial equation
      /* Define initial state of the body (rest) */
        z = 0;
        v_z = 0;
      
      equation
        z = flange_b.s;
        v_z = der(z);   /* Heave velocity */
        a_z = der(v_z); /* Heave acceleration */
        
        (a0*der(der(F_rad))) + (a1*der(F_rad)) + (a2*F_rad) = (b0*der(v_z)) + (b1*v_z); /* Transfer function representation of radiation force */
        ((M+Ainf)*a_z) + F_rad + (C*z) = F_wav; /* Cummins equation */
        
    end TF2DEModel;
     
    model SS2DEModel
      extends Modelica.Blocks.Icons.Block;
        /* State space implementation of Cummins equation to represent floating bodies
           State space representation of the radiation force convolution integral
           Using this component requires the SS Matrices to be provided in the input MATLAB struct file */
  
      /* Input interface for the wave excitation force */
      Modelica.Blocks.Interfaces.RealInput F_wav "Input flange transferring wave excitation force" annotation(
        Placement(transformation(extent = {{-90,-10},{-110,10}})));
  
      /* Output interface for the body heave dynamics and radiation force */
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b "Output flange transferring force and plant motion" annotation(
        Placement(transformation(extent = {{90,-10},{110,10}})));
  
      /* Component parameters */
      parameter Modelica.SIunits.Mass M = 727010 "Total mass of the body (including ballast)";
      parameter Modelica.SIunits.Mass Ainf = 1232838 "Added mass at maximum (cut-off) frequency";
      parameter Modelica.SIunits.TranslationalSpringConstant C = 2800951.2 "Hydrostatic stiffness";
  
      /* Component variables */
      Real A1[2,2] = {{0,1},{-1.011,-0.9366}};
      Real B1[2] = {683237,-585411};
      Real C1[2] = {1,0};
      Real D1 = 0;
      Real x1; /* Dummy variable 1 */
      Real x2; /* Dummy variable 2 */
      
      Modelica.SIunits.Length z "Vertical motion of the plant";
      Modelica.SIunits.Velocity v_z "Vertical velocity of the plant";
      Modelica.SIunits.Acceleration a_z "Vertical acceleration of the plant";
      Modelica.SIunits.Force F_rad "Radiation Force";
          
      initial equation
      /* Define initial state of the system (rest) */
        z = 0;
        v_z = 0;
      
      equation
        z = flange_b.s; /* Heave displacement */
        v_z = der(z);   /* Heave velocity */
        a_z = der(v_z); /* Heave acceleration */
        
        der(x1) = (A1[1,1]*x1) + (A1[1,2]*x2) + (B1[1]*v_z);  /* Manually define SS first matrix operation */
        der(x2) = (A1[2,1]*x1) + (A1[2,2]*x2) + (B1[2]*v_z);  /* Manually define SS second matrix operation */
        F_rad = (C1[1]*x1) + (C1[2]*x2) + (D1*v_z);           /* Assemble radiation force */
        ((M+Ainf)*a_z) + F_rad + (C*z) = F_wav;               /* Cummins equation */
      
    end SS2DEModel;
  end RadForceModels;

  model CumminsPlant
    extends Modelica.Blocks.Icons.Block;
    /* Component to solve the Cummins equation for a floating single body in heave
       State space representation of the radiation force convolution integral
       Imports hydrodynamic data and state-space matrices from an input MATLAB struct */

    /* Input interface for the wave excitation force */
    Modelica.Blocks.Interfaces.RealInput F_wav "Continuous, real wave excitation force" annotation(
      Placement(transformation(extent = {{-90,-10},{-110,10}})));

    /* Output interface for the heave dynamics and radiation force */
    Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a "Mechanical flange - displacement & force" annotation(
      Placement(transformation(extent = {{90,-10},{110,10}})));

    /* Custom connector for the heave velocity and radiation force */
    OET.Connectors.DataCollector conn() "Connector for velocity and radiation force" annotation(
        Placement(transformation(extent = {{-10,90},{10,110}})));

    /* Component parameters */
    parameter Modelica.SIunits.Mass M = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.m33",1,1)) "Total mass of the body (including ballast)";
    parameter Modelica.SIunits.Mass Ainf = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.Ainf33",1,1)) "Added mass at maximum (cut-off) frequency";
    parameter Modelica.SIunits.TranslationalSpringConstant C = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.Khs33",1,1)) "Hydrostatic stiffness";

    /* Component parameters - State space matrices import */
    parameter Real A1[2,2] = Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.ss_rad33.A",2,2) "State matrix";
    parameter Real B1[1,2] = transpose(Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.ss_rad33.B",2,1)) "Input matrix";
    parameter Real C1[1,2] = Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.ss_rad33.C",1,2) "Output matrix";
    parameter Real D1 = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.ss_rad33.D",1,1)) "Feedthrough / feedforward matrix";

    /* Component variables */
    Real x[2,1];  /* Dummy variable for state-space representation */
    Real x1;      /* Dummy variable 1 for manual state-space matrix operation */
    Real x2;      /* Dummy variable 2 for manual state-space matrix operation */
    
    Modelica.SIunits.Length z "Heave displacement";
    Modelica.SIunits.Velocity v_z "Heave velocity";
    Modelica.SIunits.Acceleration a_z "Heave acceleration";
    Modelica.SIunits.Force F_rad "Radiation Force";
        
    initial equation
    /* Define initial state of the body (rest) */
      z = 0;
      v_z = 0;
    
    equation
      flange_a.f = F_wav;  /* Wave excitation force */
      z = flange_a.s;      /* Heave displacement */
      v_z = der(z);        /* Heave velocity */
      a_z = der(v_z);      /* Heave acceleration */
      
      /* Scalable matrix state-space representation for radiation convolution */
      /* Experimental feature - Commented out */
/*      der(x) = A1*x + B1*v_z;
      F_rad = C1*x + D1*v_z; */
      
      /* Manual state-space matrix representation of radiation convolution */
      /* In lieu of above commented-out code block with a scalable matrix representation */
      der(x1) = (A1[1,1]*x1) + (A1[1,2]*x2) + (B1[1,1]*v_z);
      der(x2) = (A1[2,1]*x1) + (A1[2,2]*x2) + (B1[1,2]*v_z);
      F_rad = (C1[1,1]*x1) + (C1[1,2]*x2) + (D1*v_z);
      
      ((M+Ainf)*a_z) + F_rad + (C*z) = F_wav; /* Cummins equation */

      /* Connector definitions */
      conn.F_rad = F_rad;
      conn.v_z = v_z;
      
  end CumminsPlant;
  
  package ReferenceModels
  /* Package with reference and validation models
     Validation against WEC-Sim results
     Comparison against a frequency-dependent approach */
    
    model RefData1
    /* Validation model against WEC-Sim results imported as a MATLAB struct */

      /* Interface for heave dynamics and radiation force */
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a "Mechanical flange - displacement & force" annotation(
        Placement(transformation(extent = {{-90,-10},{-110,10}})));

      /* Custom connector for the heave velocity and radiation force */
      OET.Connectors.DataCollector conn() "Connector for velocity and radiation force" annotation(
        Placement(transformation(extent = {{-10,90},{10,110}})));

      /* Component parameters */
      parameter Real Ts = 0.1 "Discrete sampling time interval for reference dataset in seconds";
      parameter Real Tsim = 30 "Simulation duration in seconds";
      
      /* Component variables */
      Real F_rad "Radiation force";
      Real v_z "Heave velocity";
      Real RMS_error_Frad "Error in radiation force";
      Real RMS_error_v "Error in heave velocity";
      Real F_comp "Radiation force from Cummins' equation solution";
      Real v_comp "Heave velocity from Cummins' equation solution";
           
      Real F_rad_disc2[:,1] = Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.FRadRef",5001,1) "5001-sample discrete radiation force";
      Real F_rad_disc[(Tsim/Ts) + 1,1] = F_rad_disc2[1:(Tsim/Ts)+1,:] "Radiation force subset from Tsim";
      
      Real v_z_disc2[:,1] = Modelica.Utilities.Streams.readRealMatrix("C:/.../bemData.mat","bemData.vRef",5001,1) "5001-sample discrete heave velocity";
      Real v_z_disc[(Tsim/Ts) + 1,1] = v_z_disc2[1:(Tsim/Ts)+1,:] "Heave velocity subset from Tsim";
      
      equation
        /* Convert discrete validation heave radiation force to a continuous time variable */
        for k in 1:size(F_rad_disc,1) loop
          when time > Ts*(k-1)-Ts/2 and time <= Ts*k-Ts/2 then
            F_rad = F_rad_disc[k,1];
          end when;
        end for;
        
        /* Convert discrete validation heave velocity to a continuous time variable */
        for k in 1:size(v_z_disc,1) loop
          when time > Ts*(k-1)-Ts/2 and time <= Ts*k-Ts/2 then
            v_z = v_z_disc[k,1];
          end when;
        end for;
        
        ABS_error_v = abs(v_comp - v_z);      /* Absolute heave velocity error */
        ABS_error_Frad = abs(F_comp - F_rad); /* Absolute radiation force error */
        
        /* Connector declarations */
        F_comp = conn.F_rad;
        v_comp = conn.v_z;
        
    end RefData1;
    
    model FreqIndComp
    extends Modelica.Blocks.Icons.Block;
      /* Comparison model with a frequency independent representation of dynamics for a single floating body
         Modelled using the NTNU OpenModelica library (Viswanathan & Holden, 2019)
         Equation: M*der(v_z) + C/2*der(z)*(SSE_X-z3) + K*z = rho_w*g*pi*r^2*SSE_X;
         Uses value of frequency-dependent C at frequency for peak time period */

      /* Input interface for wave excitation force*/
      Modelica.Blocks.Interfaces.RealInput F_wav "Continuous, real wave excitation force signal" annotation(
      Placement(transformation(extent = {{-90,-10},{-110,10}})));

      /* Output interface for heave displacement and radiation force */
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b "Mechanical flange - displacement & force" annotation(
      Placement(transformation(extent = {{90,-10},{110,10}})));

      /* Custom connector */
      OET.Connectors.DataCollector conn() "Connector - radiation force & velocity" annotation(
        Placement(transformation(extent = {{-10,90},{10,110}})));

      /* Component parameters */
      parameter Modelica.SIunits.Mass M = 727010 "Total mass of the plant (including ballast)";
      parameter Modelica.SIunits.Mass Ainf = 1232838 "Added mass at infinite (cut-off) frequency - deprecated parameter";
      parameter Modelica.SIunits.TranslationalSpringConstant Khs = 2800951.2 "Hydrostatic restoring constant";
      parameter Modelica.SIunits.TranslationalDampingConstant C = 584473.5 "Hydrodynamic damping constant";
      constant Real g = Modelica.Constants.g_n "Value of gravitional acceleration constant";
      constant Real pi = Modelica.Constants.pi "Value of pi";
      parameter Real A = 1453672.6 "Added mass";
      parameter Real rho_w = 1025 "Water density, default for salt water @ 25 deg. Celsius";
      parameter Real T_int = 8 "Mean time period of system";

      /* Component variables */
      Modelica.SIunits.Distance z "Heave displacement";
      Modelica.SIunits.Velocity v_z "Heave velocity";
      Modelica.SIunits.Acceleration a_z "Heave acceleration";
      Modelica.SIunits.Length h_float = 3 "Submerged height of plant";
      Modelica.SIunits.Force F_rad "Radiation force";
            
      equation
        v_z = der(z);   /* Heave velocity */
        a_z = der(v_z); /* Heave acceleration */
        
        F_rad = (C*v_z*h_float);               /* Frequency-independent radiation force */
        (M+A)*a_z + F_rad + (Khs*z) = F_wav;   /* Frequency-independent representation of equation of motion */
        
        /* Connector declarations */
        conn.F_rad = F_rad;
        conn.v_z = v_z;
        
    end FreqIndComp;
  end ReferenceModels;
  
  package Functions
  /* Package defining the functions used in this library
     No functions currently - Deprecated package */
  end Functions;
  
  package Connectors
  /* Package defining custom connectors used in this library
     Using Modelica interfaces is preferred over custom-built connectors since interfacing components leverages in-built behaviour */
      
    connector DataCollector
    /* Connector for user-defined interfacing of the heave velocity and radiation force */
      Modelica.Blocks.Interfaces.RealOutput F_rad;
      Modelica.Blocks.Interfaces.RealOutput v_z;
    end DataCollector;
  end Connectors;

end OET;
/* End of library source code */

/*  Modelica ocean engineering library in MapleSim
    Developed at:
          Sys-MoDEL
          University of New Brunswick, Fredericton
          New Brunswick, E3B 5A3, Canada
    Copyright under the terms of the GNU General Public License
*/
