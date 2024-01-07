/*  Modelica Ocean Engineering Toolbox v0.2
    Copyright (C) 2024  Ajay Menon

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
     Ocean Engineering Toolbox (LIBRARY)
     |  
     |->  Wave Profile (PACKAGE)
     |    |-> Regular Wave (PACKAGE)
     |    |   |-> Regular Airy Wave (MODEL)               [!- Monochromatic wave from linear Airy wave theory]
     |    |
     |    |-> Irregular Wave (PACKAGE)
     |    |   |-> Pierson Moskowitz Spectrum (MODEL)      [!- Fully-developed sea state]
     |    |   |-> Brettschneider Spectrum (MODEL)         [!- Modified PM spectrum for developing sea state]
     |    |   |-> JONSWAP Spectrum (MODEL)                [!- Developing sea state with limited fetch]
     |    |
     |->  Structures (MODEL)
     |    |-> Plant Body (MODEL)                          [!- Solves motion through the Cummins equation]
     |
     |-> Comparison (PACKAGE)
     |    |-> Reference Dataset (MODEL)                   [!- Imported from a Matlab struct]
     |    |-> Frequency-Independent (MODEL)               [!- For comparison with the frequency-independent solution]
     |
     |->  Functions (PACKAGE)
     |    |-> waveNumber (FUNCTION)                       [!- Wave number iterations from frequency and depth]
     |    |-> randomNumberGen (FUNCTION)                  [!- Random numbers through XOR shift generator]
     |    |-> frequencySelector (FUNCTION)                [!- Select wave frequencies from a range]
     |    |-> spectrumGenerator_PM (FUNCTION)             [!- Generate Pierson Moskowitz spectrum for frequency components]
     |    |-> spectrumGenerator_BRT (FUNCTION)            [!- Generate Brettschneider spectrum for frequency components]
     |    |-> spectrumGenerator_JONSWAP (FUNCTION)        [!- Generate JONSWAP spectrum for frequency components]
     |
     |->  Connectors (PACKAGE)
          |-> WaveOutConn (CONNECTOR)                     [!- Output transfer wave elevation and excitation force]
          |-> WaveInConn (CONNECTOR)                      [!- Input transfer wave elevation and excitation force]
          |-> DataCollector (CONNECTOR)                   [!- Transfer WEC dynamics - velocity and radiation force]
*/

/* Source Code - START */
package OceanEngineeringToolbox
  /*  Library visibility */
  extends Modelica.Icons.Package;

  package WaveProfile
    package RegularWave
      model AiryWave
        /*  Airy Wave - Wave elevation profile and excitation force */
        extends Modelica.Blocks.Icons.Block;
        import Modelica.Math.Vectors;
        Modelica.Blocks.Interfaces.RealOutput F_exc "Wave excitation time series" annotation(
          Placement(transformation(extent = {{100, -10}, {120, 10}})));
        
        /*  Variable declarations */
        Real Fexc_Re[1, :] = Modelica.Utilities.Streams.readRealMatrix("F:/WEC-Sim Simulations/Float Line Attenuator/PlateBEM.mat", "PlateBEM.FexcRe", 1, 260);
        Real Fexc_Im[1, :] = Modelica.Utilities.Streams.readRealMatrix("F:/WEC-Sim Simulations/Float Line Attenuator/PlateBEM.mat", "PlateBEM.FexcIm", 1, 260);
        Real w[1, :] = Modelica.Utilities.Streams.readRealMatrix("F:/WEC-Sim Simulations/Float Line Attenuator/PlateBEM.mat", "PlateBEM.w", 1, 260);
        Real F_excRe[260] "Real component of excitation coefficient from BEMIO";
        Real F_excIm[260] "Imaginary component of excitation coefficient from BEMIO";
        Real w2[260] "Frequency distribution from BEMIO output file";
        
        constant Real pi = Modelica.Constants.pi "Mathematical constant pi";
        constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
        
        parameter Modelica.Units.SI.Length d = 100 "Water depth";
        parameter Modelica.Units.SI.Density rho = 1025 "Density of seawater";
        parameter Modelica.Units.SI.Length Hs = 2.5 "Significant Wave Height";
        parameter Modelica.Units.SI.AngularFrequency omega = 0.9423 "Peak spectral frequency";
        parameter Real Trmp = 200 "Interval for ramping up of waves during start phase";
        parameter Modelica.Units.SI.Length zeta = Hs/2 "Wave amplitude component";
        parameter Real Tp = 2*pi/omega "Wave period components";
        parameter Real k = OceanEngineeringToolbox.Functions.waveNumber(d, omega) "Wave number component";
        Real ExcCoeffRe "Real component of excitation coefficient for frequency components";
        Real ExcCoeffIm "Imaginary component of excitation coefficient for frequency components";
        Real zeta_rmp "Ramp value for surface elevation";
        Modelica.Units.SI.Length SSE "Sea surface elevation";
        Modelica.Units.SI.Length SSE_unramp "Unramped sea surface elevation";
      
      equation
        /* Convert Matlab-import matrices to vectors */
        for i in 1:260 loop
          F_excRe[i] = Fexc_Re[1, i];
          F_excIm[i] = Fexc_Im[1, i];
          w2[i] = w[1, i];
        end for;
    
        /* Define amplitude for each frequency component - ramp time */
        if time < Trmp then
          zeta_rmp = sin(pi/2*time/Trmp)*zeta;
        else
          zeta_rmp = zeta;
        end if;
          
        /* Interpolate excitation coefficients (Re & Im) for each frequency component */
        ExcCoeffRe = Modelica.Math.Vectors.interpolate(w2, F_excRe, omega);
        ExcCoeffIm = Modelica.Math.Vectors.interpolate(w2, F_excIm, omega);
        
        /* Define wave elevation profile (SSE) and excitation force */
        SSE = zeta_rmp*cos(omega*time);
        SSE_unramp = zeta*cos(omega*time);
        F_exc = ((ExcCoeffRe*zeta_rmp*cos(omega*time))-(ExcCoeffIm*zeta_rmp*sin(omega*time)))*rho*g;
        
        annotation(
          Icon(graphics = {Line(origin = {-50.91, 48.08}, points = {{-33.2809, -22.5599}, {-21.2809, -20.5599}, {-13.2809, 27.4401}, {6.71907, -20.5599}, {24.7191, -24.5599}, {42.7191, -24.5599}, {44.7191, -24.5599}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(origin = {-37, 51}, points = {{-51, 29}, {-51, -29}, {37, -29}}), Text(origin = {6, 55}, extent = {{-40, 17}, {40, -17}}, textString = "Hs"), Line(origin = {22, 4}, points = {{0, 22}, {0, -22}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-7.57, -61.12}, points = {{-82.4341, -12.8774}, {-76.4341, -2.87735}, {-72.4341, -6.87735}, {-62.4341, 13.1226}, {-50.4341, -26.8774}, {-46.4341, -20.8774}, {-38.4341, -26.8774}, {-34.4341, -18.8774}, {-34.4341, 3.12265}, {-26.4341, 1.12265}, {-20.4341, 7.12265}, {-12.4341, 9.12265}, {-8.43408, 19.1226}, {1.56592, -4.87735}, {7.56592, -24.8774}, {19.5659, -6.87735}, {21.5659, 9.12265}, {31.5659, 13.1226}, {39.5659, -0.87735}, {43.5659, 11.1226}, {55.5659, 15.1226}, {63.5659, 27.1226}, {79.5659, -22.8774}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Rectangle(origin = {100, 0}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.05));
      end AiryWave;
    
    end RegularWave;
  end WaveProfile;

  package Functions
    /* Package defining explicit library functions */

    function waveNumber
    /*  Function to iteratively compute the wave number from the frequency components */
      input Real d "Water depth";
      input Real omega[:] "Wave frequency components";
      output Real k[size(omega, 1)] "Wave number components";
    protected
      constant Real g = Modelica.Constants.g_n;
      constant Real pi = Modelica.Constants.pi;
      parameter Integer n = size(omega, 1);
      Real T[size(omega, 1)] "Wave period components";
      Real L0[size(omega, 1)] "Deepwater wave length";
      Real L1(start = 0, fixed = true) "Temporary variable";
      Real L1c(start = 0, fixed = true) "Temporary variable";
      Real L[size(omega, 1)] "Iterated wave length";
    
    algorithm
      T := 2*pi./omega;
      L0 := g*T.^2/(2*pi);
      for i in 1:size(omega, 1) loop
        L1 := L0[i];
        L1c := 0;
        while abs(L1c - L1) > 0.001 loop
          L1c := L1;
          L[i] := g*T[i]^2/(2*pi)*tanh(2*pi/L1*d);
          L1 := L[i];
        end while;
      end for;
      k := 2*pi./L;
    end waveNumber;

  end Functions;

  package Connectors
  /*  Package defining library connectors between models */
    
    connector WaveOutConn
    /*  Output datastream - wave elevation & excitation force */
      //Modelica.Blocks.Interfaces.RealOutput SSE;
      Modelica.Blocks.Interfaces.RealOutput F_exc;
    end WaveOutConn;
    
    connector WaveInConn
    /*  Input datastream - wave elevation & excitation force */
      //Modelica.Blocks.Interfaces.RealInput SSE;
      Modelica.Blocks.Interfaces.RealInput F_exc;
    end WaveInConn;
    
    connector DataCollector
    /*  Output datastream - velocity and radiation force */
      Modelica.Blocks.Interfaces.RealOutput F_rad;
      Modelica.Blocks.Interfaces.RealOutput v_z;
    end DataCollector;
    
  end Connectors;

end OceanEngineeringToolbox;
/* Source Code - END */

/*  Modelica Ocean Engineering Toolbox (OET)
    Developed at:
          Sys-MoDEL, 
          University of New Brunswick, Fredericton
          New Brunswick, E3B 5A3, Canada
    Copyright under the terms of the GNU General Public License
*/
