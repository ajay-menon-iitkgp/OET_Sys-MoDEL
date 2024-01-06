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
   
   WECHydrodynamics (LIBRARY)
   |->  Excitation Force (PACKAGE)
   |    |-> Regular Waves (PACKAGE)
   |    |   |-> Regular Airy Waves (MODEL) [!- Mono- and polychromatic]
   |    |
   |    |-> Irregular Waves (PACKAGE) [!- Consider combining all these models to a single component with a drop-down selection]
   |    |   |-> Pierson Mosowitz Spectrum (MODEL)
   |    |   |-> JONSWAP Spectrum (MODEL)
   |    |   |-> Bretschneider Spectrum (MODEL)
   |    |
   |    |-> Assemble Excitation Force (MODEL)
   |
   |->  Plant Model (MODEL)
   |    |-> TF2DE (MODEL) [!- If the user has the transfer function ODE representation]
   |    |-> SS2DE (MODEL) [!- Default method of solving Cummins' equation]
   |
   |->  Reference Models (PACKAGE)
   |    |-> Ref Data 1 (MODEL) [!- Imported from Matlab struct]
   |    |-> NTNU Equivalent [!- Only for comparison purposes]
   |
   |->  Functions (PACKAGE) - DEPRECATED
   |    |-> discreteToContinuous (FUNCTION)
   |
   |->  Connectors (PACKAGE)
        |-> Data Collector (CONNECTOR) [!- Transfer non-flange data]
*/

/* Important notes:
   
   1. To change the library name, replace all instances of "WECHydrodynamics" with your chosen name.
   2. The statement "extends Modelica.Blocks.Icons.Block;" is used to define visibility of the model as a block in the Modelica solver. Removing this statement causes the relevant model to disappear from the library.
   3. This library can be run in any Modelica solver engine, although it is intended for use in MapleSim. Annotations have been prepared accordingly.

    Present capabilities:
    1. The present library is capable of importing data from a struct in MATLAB to define the hydrodynamic coefficients, excitation force, and state-space quadruple
    2. Capable of estimating the displacement, velocity, acceleration, and radiation forces acting on a single body defined using the state space quadruple
*/

/*  How to use the library:
    1.  Download this program. It can be saved at any location, and does not need to specifically be in MapleSim's or Modelica's main directory.
    2.  Run a BEM solver and obtain the hydrodynamic coefficients and excitation force as a discrete array against simulation time.
    3.  The list of present requirements from the BEM solver are:
          (i)   Mass (m33)
          (ii)  Infinite added mass (Ainf33)
          (iii) Hydrostatic restoring coefficient (Khs33)
          (iv)  Wave excitation force (F_wav)
          (v)   Reference velocity from validation dataset (vRef)
          (vi)  Reference radiation force from validation dataset (FRadRef)
          (vii) State Space quadruple:
                    (a) A  (b) B  (c) C  (d) D
    4.  The results of the BEM solver must be stored in a single MATLAB struct and saved to any location.
    5.  In the Modelica code, identify the code "Modelica.Utilities.Streams" and replace the address with the location of the BEM results .mat file in the specified format. Note, the complete address must be entered, and not merely the file's name.
    6.  Import the library into MapleSim using "File -> Modelica -> Import library and documentation".
    7.  Connect the ExcitationForce model to the left interface of the Plant model.
    8.  If you wish to validate the solver's results, connect the right Flange interface of SS2DE to the RefData model.
    9.  Set the simulation time in the MapleSim toolbar, as well as on the parameters panel of ExctitationForce and RefData (Tsim).
    10. Attach probes to monitor the required variables and run the simulation. To view radiation force, attach a probe to the component and NOT the flange. From the flange, only mechanical properties can be observed- displacement (s), velocity (v), acceleration (a), force on the flange (f).
*/

// Source code
// Library name: WEC Hydrodynamics
package WECHydrodynamics

  package ExcitationForce
  // Package defning the excitation force acting on the plant model
  // The excitation force is defined as the convolution of the excitation coefficient Tau with the incident wave spectrum
  // For trial purposes, a defined excitation force is also included
  // Users can opt for regular monochromatic/polychromatic or irregular spectrums
  // The irregular waves package currently supports the Pierson Moskowitz spectrum, but can include JONSWAP & Bretschneider spectra as well
  
    package RegularWaves
    // Package defining regular wave spectrums - UNDER CONSTRUCTION
    // Defined for both monochromatic & polychromatic spectra
    // Users control parameters of the wave spectrum
    // The components of this package output the surface elevation
      extends Modelica.Blocks.Icons.Block;
      
      Modelica.Blocks.Interfaces.RealOutput eta "Continuous, real excitation force signal"
      annotation(
        Placement(
          transformation(
            extent = {
              {100,-10},
              {120,10}
            }
          )
        )
      );
      
      parameter Real n_components = 1;
      parameter Modelica.SIunits.AngularFrequency w_min = 0.5;
      parameter Modelica.SIunits.AngularFrequency w_max = 1.1;
      
      Modelica.SIunits.AngularFrequency w;
      Modelica.SIunits.Length A;
      Modelica.SIunits.WaveNumber k;
      equation
        
        
      
    end RegularWaves;
    
    package IrregularWaves
    // Package defining irregular wave spectrums - UNDER CONSTRUCTION
    // Defined for the Pierson Moskowitz spectrum, can be expanded to the Bretschneider & JONSWAP spectra as well
    // Users define the irregular spectral parameters
    // The components of this package output the surface elevation
      
    end IrregularWaves;
    
    package UserDefinedWave
    // Package defining wave elevation from user input - UNDER CONSTRUCTION
    // Defined for any random discrete wave elevation/sea spectrum
    // Users specify the location of the MATLAB file storing the discrete wave elevation distribution
    // The components of this package output the continuous surface elevation
      
    end UserDefinedWave;
  
    model AssembleExcitationForce
    // Component that generates the excitation force from the wave elevation
    // Inputs the wave surface elevations from the regular/irregular/user-defined packages
    // Performs convolution with the excitation impulse response function
    // Requires the excitation convolution state-space matrices to be specified
    // Yields the excitation force for the wave elevation
      extends Modelica.Blocks.Icons.Block;
      
      Modelica.Blocks.Interfaces.RealInput eta "Wave elevation profile"
      annotation(
        Placement(
          transformation(
            extent = {
              {-100,-10},
              {-120,10}
            }
          )
        )
      );
      
      Modelica.Blocks.Interfaces.RealOutput F_exc "Continuous, real excitation force signal"
      annotation(
        Placement(
          transformation(
            extent = {
              {100,-10},
              {120,10}
            }
          )
        )
      );
      
      parameter Real A1[2,2] = Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.ss_exc.A",2,2) "State matrix";
      parameter Real B1[1,2] = transpose(Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.ss_exc.B",2,1)) "Input matrix";
      parameter Real C1[1,2] = Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.ss_exc.C",1,2) "Output matrix";
      parameter Real D1 = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.ss_exc.D",1,1)) "Feedthrough / feedforward matrix";
      
      Real x1;
      Real x2;
      
      equation
        der(x1) = (A1[1,1]*x1) + (A1[1,2]*x2) + (B1[1,1]*eta);
        der(x2) = (A1[2,1]*x1) + (A1[2,2]*x2) + (B1[1,2]*eta);
        F_exc = (C1[1,1]*x1) + (C1[1,2]*x2) + (D1*eta);
        
    end AssembleExcitationForce;
    
    model ImportExcitationForce
    // Model to import the sample excitation force from a MATLAB struct
    // Uses the Utilities.Streams library to import MATLAB matrix data
    // Location of the MATLAB struct file must be specified in the code
    // Converts the discrete import samples to a continuous excitation force
      extends Modelica.Blocks.Icons.Block;
      
      Modelica.Blocks.Interfaces.RealOutput F_exc2 "Continuous, real excitation force signal"
      annotation(
        Placement(
          transformation(
            extent = {
              {100,-10},
              {120,10}
            }
          )
        )
      );
      
      parameter Modelica.SIunits.Time Ts = 0.1 "Discrete sampling time interval in seconds";
      parameter Modelica.SIunits.Time Tsim = 30 "Simulation duration in seconds";
      
      Real F_exc_discrete2[:,1] = Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.F_wav",5001,1) "5001-samples discrete excitation force";
      Real F_exc_discrete[(Tsim/Ts) + 1,1] = F_exc_discrete2[1:(Tsim/Ts)+1,:] "Discret excitation force subset based on desired simulation duration";
      
      equation
        // Loop through the discrete subset
        // When time is between half time-steps, hold the discrete value
        // Builds the continuous function as a step-function of the discrete samples
        for k in 1:size(F_exc_discrete,1) loop
          when time > Ts*(k-1)-Ts/2 and time <= Ts*k-Ts/2 then
            F_exc2 = F_exc_discrete[k,1];
          end when;
        end for;
      
    end ImportExcitationForce;
     
  end ExcitationForce;
    
  model TF2DEModel
  // Transfer function to differential equation approach
  // Models the plant using Cummins' equation
  // Uses the transfer function to solve for the radiation convolution
  // Requires the TF's differential equation to be coded-in
    extends Modelica.Blocks.Icons.Block;
    
    Modelica.Blocks.Interfaces.RealInput F_wav "Input flange transferring wave excitation force"
    annotation(
      Placement(
        transformation(
          extent = {
            {-90,-10},
            {-110,10}
          }
        )
      )
    );
    
    Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b "Output flange transferring force and plant motion"
    annotation(
      Placement(
        transformation(
          extent = {
            {90,-10},
            {110,10}
          }
        )
      )
    );
      
    parameter Modelica.SIunits.Mass M = 727010 "Total mass of the body (including ballast)";
    parameter Modelica.SIunits.Mass Ainf = 1232838 "Added mass at maximum (cut-off) frequency";
    parameter Modelica.SIunits.TranslationalSpringConstant C = 2800951.2 "Hydrostatic stiffness";
    parameter Real a0 = 1 "Coefficient a0 of second-order radiation force derivative";
    parameter Real a1 = 0.9366 "Coefficient a1 of first-order radiation force derivative";
    parameter Real a2 = 1.011 "Coefficient a2 of zeroth-order radiation force derivative";
    parameter Real b0 = 683237 "Coefficient b0 of first-order velocity derivative";
    parameter Real b1 = 54478 "Coefficient b1 of zeroth-order velocity derivative";
    Modelica.SIunits.Length z "Heave displacement";
    Modelica.SIunits.Velocity v_z "Heave velocity";
    Modelica.SIunits.Acceleration a_z "Heave acceleration";
    Modelica.SIunits.Force F_rad "Radiation Force";
    
    initial equation
    // Define body at initial rest
      z = 0;
      v_z = 0;
    
    equation
      z = flange_b.s; // Plant motion is tied to the vertical displacement
//      flange_b.f = F_wav; // Force acting on the plant is excitation force
      v_z = der(z);   // Heave velocity
      a_z = der(v_z); // Heave acceleration
      
      (a0*der(der(F_rad))) + (a1*der(F_rad)) + (a2*F_rad) = (b0*der(v_z)) + (b1*v_z); // Transfer function representation of radiation force
      ((M+Ainf)*a_z) + F_rad + (C*z) = F_wav; // Cummins' equation
      
  end TF2DEModel;
   
  model SS2DEModel
  // Plant model to solve Cummins' equation
  // Represent radiation convolution using the state-space approach
  // State-space quadruple represented as ODE
  // Matrix approach allows for scalability of the system
    extends Modelica.Blocks.Icons.Block;
    
    Modelica.Blocks.Interfaces.RealInput F_wav "Input flange transferring wave excitation force"
    annotation(
      Placement(
        transformation(
          extent = {
            {-90,-10},
            {-110,10}
          }
        )
      )
    );
    
    Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b "Output flange transferring force and plant motion"
    annotation(
      Placement(
        transformation(
          extent = {
            {90,-10},
            {110,10}
          }
        )
      )
    );
      
    parameter Modelica.SIunits.Mass M = 727010 "Total mass of the body (including ballast)";
    parameter Modelica.SIunits.Mass Ainf = 1232838 "Added mass at maximum (cut-off) frequency";
    parameter Modelica.SIunits.TranslationalSpringConstant C = 2800951.2 "Hydrostatic stiffness";
    
    Real A1[2,2] = {{0,1},{-1.011,-0.9366}};
    Real B1[2] = {683237,-585411};
    Real C1[2] = {1,0};
    Real D1 = 0;
    Real x1;  // Dummy 1
    Real x2;  // Dummy 2
    
    Modelica.SIunits.Length z "Vertical motion of the plant";
    Modelica.SIunits.Velocity v_z "Vertical velocity of the plant";
    Modelica.SIunits.Acceleration a_z "Vertical acceleration of the plant";
    Modelica.SIunits.Force F_rad "Radiation Force";
        
    initial equation
      z = 0;
      v_z = 0;
    
    equation
      z = flange_b.s; // Plant motion is tied to the vertical displacement
      v_z = der(z);   // Define equation for vertical plant velocity
      a_z = der(v_z); // Define equation for vertical plant acceleration
      
      der(x1) = (A1[1,1]*x1) + (A1[1,2]*x2) + (B1[1]*v_z);
      der(x2) = (A1[2,1]*x1) + (A1[2,2]*x2) + (B1[2]*v_z);
      F_rad = (C1[1]*x1) + (C1[2]*x2) + (D1*v_z);
      ((M+Ainf)*a_z) + F_rad + (C*z) = F_wav;
    
  end SS2DEModel;
  
  model CumminsPlant
  // Component to solve Cummins' equation
  // Uses the state-space approach to solve radiation convolution
  // Imports hydrodynamic data and state-space quadruple from MATLAB struct
  // Transfers velocity and radiation force through conn() connector
    extends Modelica.Blocks.Icons.Block;
    
    Modelica.Blocks.Interfaces.RealInput F_wav "Continuous, real wave excitation force"
    annotation(
      Placement(
        transformation(
          extent = {
            {-90,-10},
            {-110,10}
          }
        )
      )
    );
    
    Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a "Mechanical flange - displacement & force"
    annotation(
      Placement(
        transformation(
          extent = {
            {90,-10},
            {110,10}
          }
        )
      )
    );
    
    WECHydrodynamics.Connectors.DataCollector conn() "Connector for velocity and radiation force"
      annotation(
        Placement(
          transformation(
            extent = {
              {-10,90},
              {10,110}
            }
          )
        )
      );
      
    parameter Modelica.SIunits.Mass M = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.m33",1,1)) "Total mass of the body (including ballast)";
    parameter Modelica.SIunits.Mass Ainf = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.Ainf33",1,1)) "Added mass at maximum (cut-off) frequency";
    parameter Modelica.SIunits.TranslationalSpringConstant C = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.Khs33",1,1)) "Hydrostatic stiffness";
    
    parameter Real A1[2,2] = Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.ss_rad33.A",2,2) "State matrix";
    parameter Real B1[1,2] = transpose(Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.ss_rad33.B",2,1)) "Input matrix";
    parameter Real C1[1,2] = Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.ss_rad33.C",1,2) "Output matrix";
    parameter Real D1 = scalar(Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.ss_rad33.D",1,1)) "Feedthrough / feedforward matrix";
    
    Real x[2,1]; // Dummy variables for state-space representation
    Real x1;
    Real x2;
    
    Modelica.SIunits.Length z "Heave displacement";
    Modelica.SIunits.Velocity v_z "Heave velocity";
    Modelica.SIunits.Acceleration a_z "Heave acceleration";
    Modelica.SIunits.Force F_rad "Radiation Force";
        
    initial equation
    // Define body at rest initially
      z = 0;
      v_z = 0;
    
    equation
      flange_a.f = F_wav; // Excitation force acts on plant
      z = flange_a.s; // Plant motion is tied to the vertical displacement
      v_z = der(z);   // Heave velocity
      a_z = der(v_z); // Heave acceleration
      
      // Scalable matrix of state-space representation for radiation convolution
      // Deprecated
/*      der(x) = A1*x + B1*v_z;
      F_rad = C1*x + D1*v_z;*/
      
      // Broken-down state-space representation of radiation convolution
      // If above code returns error or incorrect results, comment the above section out and replace with below comment block
      der(x1) = (A1[1,1]*x1) + (A1[1,2]*x2) + (B1[1,1]*v_z);
      der(x2) = (A1[2,1]*x1) + (A1[2,2]*x2) + (B1[1,2]*v_z);
      F_rad = (C1[1,1]*x1) + (C1[1,2]*x2) + (D1*v_z);
      
      // Cummins' equation
      ((M+Ainf)*a_z) + F_rad + (C*z) = F_wav;
      
      // Connector declarations
      conn.F_rad = F_rad;
      conn.v_z = v_z;
      
  end CumminsPlant;
    
  package ReferenceModels
  // Package defining reference models
  // For validation against WEC-Sim results
  // Includes an equivalent representation of NTNU's ocean engineering library in OpenModelica
    
    model RefData1
    // Reference model for validation against an external dataset
    // Reference imported from a MATLAB struct
    // Can be used to compare against WEC-Sim results
    
      Modelica.Mechanics.Translational.Interfaces.Flange_a flange_a "Mechanical flange - displacement & force"
      annotation(
        Placement(
          transformation(
            extent = {
              {-90,-10},
              {-110,10}
            }
          )
        )
      );
      
      WECHydrodynamics.Connectors.DataCollector conn() "Connector for velocity and radiation force"
      annotation(
        Placement(
          transformation(
            extent = {
              {-10,90},
              {10,110}
            }
          )
        )
      );
      
      Real F_rad "Radiation force";
      Real v_z "Heave velocity";
      Real RMS_error_Frad "Error in radiation force";
      Real RMS_error_v "Error in heave velocity";
      Real F_comp "Radiation force from Cummins' equation solution";
      Real v_comp "Heave velocity from Cummins' equation solution";
      
      parameter Real Ts = 0.1 "Discrete sampling time interval for reference dataset in seconds";
      parameter Real Tsim = 30 "Simulation duration in seconds";
      
      Real F_rad_disc2[:,1] = Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.FRadRef",5001,1) "5001-sample discrete radiation force";
      Real F_rad_disc[(Tsim/Ts) + 1,1] = F_rad_disc2[1:(Tsim/Ts)+1,:] "Radiation force subset from Tsim";
      
      Real v_z_disc2[:,1] = Modelica.Utilities.Streams.readRealMatrix("C:/Users/r75zs/Desktop/Matlab files/Convo 2/bemData.mat","bemData.vRef",5001,1) "5001-sample discrete heave velocity";
      Real v_z_disc[(Tsim/Ts) + 1,1] = v_z_disc2[1:(Tsim/Ts)+1,:] "Heave velocity subset from Tsim";
      
      equation
        // Convert discrete reference radiation force to continuous variable
        for k in 1:size(F_rad_disc,1) loop
          when time > Ts*(k-1)-Ts/2 and time <= Ts*k-Ts/2 then
            F_rad = F_rad_disc[k,1];
          end when;
        end for;
        
        // Convert discrete reference heave velocity to continuous variable
        for k in 1:size(v_z_disc,1) loop
          when time > Ts*(k-1)-Ts/2 and time <= Ts*k-Ts/2 then
            v_z = v_z_disc[k,1];
          end when;
        end for;
        
        // Unable to compute RMS error between the reference and computed data
        // Deprecated code
        // dummy1 = der(RMS_error_Frad^2);
        // dummy2 = der(RMS_error_v^2);
        // Tsim*dummy1 = ((F_comp - F_rad)^2);
        // Tsim*dummy2 = ((v_comp - v_z)^2);
        RMS_error_v = abs(v_comp - v_z);      // Absolute heave velocity error
        RMS_error_Frad = abs(F_comp - F_rad); // Absolute radiation force error
        
        // Connector declarations
        F_comp = conn.F_rad;
        v_comp = conn.v_z;
        
    end RefData1;
    
    model NTNURefData
    // Reference model equivalent of NTNU's OpenModelica library
    // Uses similar concepts: M*der(v_z) + C/2*der(z)*(SSE_X-z3) + K*z = rho_w*g*pi*r^2*SSE_X;
    // Using value of frequency-dependent C at frequency for mean time-period
    // Intended to be used to compare improvement of the WECHydrodynamics library over NTNU's library
      extends Modelica.Blocks.Icons.Block;
      
      Modelica.Blocks.Interfaces.RealInput F_wav "Continuous, real wave excitation force signal"
      annotation(
      Placement(
        transformation(
          extent = {
            {-90,-10},
            {-110,10}
          }
        )
      )
    );
    
      Modelica.Mechanics.Translational.Interfaces.Flange_b flange_b "Mechanical flange - displacement & force"
      annotation(
      Placement(
        transformation(
          extent = {
            {90,-10},
            {110,10}
          }
        )
      )
    );
      
      WECHydrodynamics.Connectors.DataCollector conn() "Connector - radiation force & velocity"
      annotation(
        Placement(
          transformation(
            extent = {
              {-10,90},
              {10,110}
            }
          )
        )
      );
      
      parameter Modelica.SIunits.Mass M = 727010 "Total mass of the plant (including ballast)";
      parameter Modelica.SIunits.Mass Ainf = 1232838 "Added mass at infinite (cut-off) frequency - deprecated parameter";
      parameter Modelica.SIunits.TranslationalSpringConstant Khs = 2800951.2 "Hydrostatic restoring constant";
      parameter Modelica.SIunits.TranslationalDampingConstant C = 584473.5 "Hydrodynamic damping constant";
      constant Real g = Modelica.Constants.g_n "Value of gravitional acceleration constant";
      constant Real pi = Modelica.Constants.pi "Value of pi";
      
      parameter Real A = 1453672.6 "Added mass";
      parameter Real rho_w = 1025 "Water density, default for salt water @ 25 deg. Celsius";
      parameter Real T_int = 8 "Mean time period of system";
      
      Modelica.SIunits.Distance z "Heave displacement";
      Modelica.SIunits.Velocity v_z "Heave velocity";
      Modelica.SIunits.Acceleration a_z "Heave acceleration";
      
      Modelica.SIunits.Length h_float = 3 "Submerged height of plant";
      Modelica.SIunits.Force F_rad "Radiation force";
            
      equation
        v_z = der(z); // Heave velocity
        a_z = der(v_z); // Heave acceleration
        
        F_rad = (C*v_z*h_float); // Linear radiation force
        (M+A)*a_z + F_rad + (Khs*z) = F_wav;  // Linearized representation of plant's equation of motion
        
        // Connector declaration
        conn.F_rad = F_rad;
        conn.v_z = v_z;
        
    end NTNURefData;
        
  end ReferenceModels;
  
  package Functions
  // Package defining the functions used in this library
  // None at the moment - Deprecated package
    
    function ConvIntegral2
    // Deprecated functin
      input Real K[:];
      input Real Ts;
      input Real u[:];
      output Integer F_rad;

      algorithm
        F_rad := sum(u .* K) * Ts;
        
    end ConvIntegral2;
    
  end Functions;
  
  package Connectors
  // Package defining the connectors used in this library
  // Preference is to use Modelica interfaces over custom-built Connectors
  // Interfacing components leverages in-built behaviour and is convenient to use
  // Force components (excitation & radiation) are transferred to the plant through flanges
  // Plant motion is solved from the plant model and represented from the position property of the flange interface
      
    connector DataCollector
    // Connector for user-defined interfacing
    // Velocity and radiation force in datastream
      Modelica.Blocks.Interfaces.RealOutput F_rad;
      Modelica.Blocks.Interfaces.RealOutput v_z;
    end DataCollector;
        
  end Connectors;

end WECHydrodynamics;
// End of source code

/*  Modelica ocean engineering library in MapleSim
    Developed at:
          Sys-MoDEL
          University of New Brunswick, Fredericton
          New Brunswick, E3B 5A3, Canada
    Copyright under the terms of the GNU General Public License
*/
