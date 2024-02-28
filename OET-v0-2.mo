/*  Modelica Ocean Engineering Toolbox v0.2
    Copyright (C) 2024  Ajay Menon, Ali Haider, Kush Bubbar

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
  
/*  Library structure:
     Ocean Engineering Toolbox (LIBRARY)
     |  
     |->  Wave Profile (PACKAGE)
     |    |-> Regular Wave (PACKAGE)
     |    |   |-> LinearModel (MODEL)               [!- Monochromatic regular wave]
     |    |
     |    |-> Irregular Wave (PACKAGE)
     |    |   |-> Pierson Moskowitz Spectrum (MODEL)      [!- Fully-developed sea state]
     |    |   |-> Brettschneider Spectrum (MODEL)         [!- Modified PM spectrum for developing sea state]
     |    |   |-> JONSWAP Spectrum (MODEL)                [!- Developing sea state with limited fetch]
     |    |
     |->  Structures (PACKAGE)
     |    |-> WEC (MODEL)                       [!- Solves motion through the Cummins equation]
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

/* Source Code */

package OceanEngineeringToolbox
  extends Modelica.Icons.Package;   /*  Library visibility */

  package WaveProfile
  
    package RegularWave
      /* Package for regular wave elevation profile and excitation force calculations */

      model LinearWave
        /*  Wave elevation profile and excitation force */
        
        extends Modelica.Blocks.Icons.Block;
        import Modelica.Math.Vectors;
        OceanEngineeringToolbox.Internal.Connectors.WaveOutConn wconn;
        
        /* Environmental constants */
        constant Real pi = Modelica.Constants.pi "Mathematical constant pi";
        constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
        
        /*  Variable declarations */
        parameter String fileName;
        
        parameter Real F_excRe[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.FexcRe", 1, 260));
        parameter Real F_excIm[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.FexcIm", 1, 260));
        parameter Real w[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.w", 1, 260));
        
        parameter Modelica.Units.SI.Length d = 100 "Water depth";
        parameter Modelica.Units.SI.Density rho = 1025 "Density of seawater";
        parameter Modelica.Units.SI.Length Hs = 2.5 "Significant Wave Height";
        parameter Modelica.Units.SI.AngularFrequency omega = 0.9423 "Wave frequency";
        parameter Real Trmp = 100 "Interval for ramping up of waves during start phase";
        parameter Modelica.Units.SI.Length zeta = Hs/2 "Wave amplitude";
        parameter Real Tp = 2*pi/omega "Wave period";
        parameter Real k = 2*pi/(1.56*(Tp^2)) "Wave number";
        Real ExcCoeffRe "Real component of excitation coefficient";
        Real ExcCoeffIm "Imaginary component of excitation coefficient";
        Real zeta_rmp "Ramp value for surface elevation";
        Modelica.Units.SI.Length SSE "Sea surface elevation";
        Modelica.Units.SI.Length SSE_unramp "Unramped sea surface elevation";
      
      equation
        /* Define amplitude for each frequency component - ramp time */
        if time < Trmp then
          zeta_rmp = sin(pi/2*time/Trmp)*zeta;
        else
          zeta_rmp = zeta;
        end if;
        
        /* Interpolate excitation coefficients (Re & Im) for each frequency component */
        ExcCoeffRe = Modelica.Math.Vectors.interpolate(w, F_excRe, omega);
        ExcCoeffIm = Modelica.Math.Vectors.interpolate(w, F_excIm, omega);

        /* Define wave elevation profile (SSE) and excitation force */
        SSE = zeta_rmp*cos(omega*time);
        SSE_unramp = zeta*cos(omega*time);
        wconn.F_exc = ((ExcCoeffRe*zeta_rmp*cos(omega*time)) - (ExcCoeffIm*zeta_rmp*sin(omega*time)))*rho*g;
        
        annotation(
          Icon(graphics = {Line(origin = {-50.91, 48.08}, points = {{-33.2809, -22.5599}, {-21.2809, -20.5599}, {-13.2809, 27.4401}, {6.71907, -20.5599}, {24.7191, -24.5599}, {42.7191, -24.5599}, {44.7191, -24.5599}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(origin = {-37, 51}, points = {{-51, 29}, {-51, -29}, {37, -29}}), Text(origin = {6, 55}, extent = {{-40, 17}, {40, -17}}, textString = "Hs"), Line(origin = {22, 4}, points = {{0, 22}, {0, -22}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-7.57, -61.12}, points = {{-82.4341, -12.8774}, {-76.4341, -2.87735}, {-72.4341, -6.87735}, {-62.4341, 13.1226}, {-50.4341, -26.8774}, {-46.4341, -20.8774}, {-38.4341, -26.8774}, {-34.4341, -18.8774}, {-34.4341, 3.12265}, {-26.4341, 1.12265}, {-20.4341, 7.12265}, {-12.4341, 9.12265}, {-8.43408, 19.1226}, {1.56592, -4.87735}, {7.56592, -24.8774}, {19.5659, -6.87735}, {21.5659, 9.12265}, {31.5659, 13.1226}, {39.5659, -0.87735}, {43.5659, 11.1226}, {55.5659, 15.1226}, {63.5659, 27.1226}, {79.5659, -22.8774}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Rectangle(origin = {100, 0}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.05));
      end LinearWave;
    
    end RegularWave;

    package IrregularWave
      /* Package for irregular wave elevation profile and excitation force calculations */

      model PiersonMoskowitzWave
        /*  Pierson Moskowitz Spectrum - Wave elevation profile and excitation force */
        
        extends Modelica.Blocks.Icons.Block;
        import Modelica.Math.Vectors;
        OceanEngineeringToolbox.Internal.Connectors.WaveOutConn wconn;
        
        /* Modelica.Blocks.Interfaces.RealOutput F_exc "Wave time series" annotation(
                  Placement(transformation(extent = {{100, -10}, {120, 10}}))); */
        
        /* Environmental constants */
        constant Real pi = Modelica.Constants.pi "Mathematical constant pi";
        constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
        
        /*  Variable declarations */
        parameter String fileName;
        
        // Integer varSize[:] = Modelica.Utilities.Streams.readMatrixSize(fileName, "hydroCoeff.FexcRe");
        
        parameter Real F_excRe[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.FexcRe", 1, 260));
        parameter Real F_excIm[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.FexcIm", 1, 260));
        parameter Real w[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.w", 1, 260));
        
        parameter Modelica.Units.SI.Length d = 100 "Water depth";
        parameter Modelica.Units.SI.Density rho = 1025 "Density of seawater";
        parameter Modelica.Units.SI.Length Hs = 2.5 "Significant Wave Height";
        parameter Modelica.Units.SI.AngularFrequency omega_min = 0.03141 "Lowest frequency component/frequency interval";
        parameter Modelica.Units.SI.AngularFrequency omega_max = 3.141 "Highest frequency component";
        parameter Modelica.Units.SI.AngularFrequency omega_peak = 0.9423 "Peak spectral frequency";
        parameter Real spectralWidth_min = 0.07 "Lower spectral bound for JONSWAP";
        parameter Real spectralWidth_max = 0.09 "Upper spectral bound for JONSWAP";
        parameter Integer n_omega = 100 "Number of frequency components";
        parameter Integer localSeed = 614657 "Local random seed";
        parameter Integer globalSeed = 30020 "Global random seed";
        parameter Real rnd_shft[n_omega] = OceanEngineeringToolbox.Internal.Functions.randomNumberGen(localSeed, globalSeed, n_omega);
        parameter Integer localSeed1 = 614757 "Local rand seed";
        parameter Integer globalSeed1 = 40020 "Global rand seed";
        parameter Real epsilon[n_omega] = OceanEngineeringToolbox.Internal.Functions.randomNumberGen(localSeed1, globalSeed1, n_omega) "Wave components phase shift";
        parameter Real Trmp = 100 "Interval for ramping up of waves during start phase";
        parameter Real omega[n_omega] = OceanEngineeringToolbox.Internal.Functions.frequencySelector(omega_min, omega_max, rnd_shft);
        parameter Real S[n_omega] = OceanEngineeringToolbox.Internal.Functions.spectrumGenerator_PM(Hs, omega) "Spectral values for frequency components";
        parameter Modelica.Units.SI.Length zeta[n_omega] = sqrt(2*S*omega_min) "Wave amplitude component";
        parameter Real Tp[n_omega] = 2*pi./omega "Wave period components";
        parameter Real k[n_omega] = OceanEngineeringToolbox.Internal.Functions.waveNumber(d, omega) "Wave number component";
        Real ExcCoeffRe[n_omega] "Real component of excitation coefficient for frequency components";
        Real ExcCoeffIm[n_omega] "Imaginary component of excitation coefficient for frequency components";
        Real zeta_rmp[n_omega] "Ramp value for surface elevation";
        Modelica.Units.SI.Length SSE "Sea surface elevation";
        Modelica.Units.SI.Length SSE_unramp "Unramped sea surface elevation";
        
      equation

        /* Define amplitude for each frequency component - ramp time */
        for i in 1:n_omega loop
          if time < Trmp then
            zeta_rmp[i] = (1 + cos(pi + (pi*time/Trmp)))*zeta[i];
          else
            zeta_rmp[i] = zeta[i];
          end if;
          /* Interpolate excitation coefficients (Re & Im) for each frequency component */
          ExcCoeffRe[i] = Modelica.Math.Vectors.interpolate(w, F_excRe, omega[i])*rho*g;
          ExcCoeffIm[i] = Modelica.Math.Vectors.interpolate(w, F_excIm, omega[i])*rho*g;
        end for;
        
        /* Define wave elevation profile (SSE) and excitation force */
        SSE = sum(zeta_rmp.*cos(omega*time - 2*pi*epsilon));
        SSE_unramp = sum(zeta.*cos(omega*time - 2*pi*epsilon));
        wconn.F_exc = sum((ExcCoeffRe.*zeta_rmp.*cos(omega*time - 2*pi*epsilon)) - (ExcCoeffIm.*zeta_rmp.*sin(omega*time - 2*pi*epsilon)));
        
        annotation(
          Icon(graphics = {Line(origin = {-50.91, 48.08}, points = {{-33.2809, -22.5599}, {-21.2809, -20.5599}, {-13.2809, 27.4401}, {6.71907, -20.5599}, {24.7191, -24.5599}, {42.7191, -24.5599}, {44.7191, -24.5599}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(origin = {-37, 51}, points = {{-51, 29}, {-51, -29}, {37, -29}}), Text(origin = {6, 55}, extent = {{-40, 17}, {40, -17}}, textString = "Hs"), Line(origin = {22, 4}, points = {{0, 22}, {0, -22}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-7.57, -61.12}, points = {{-82.4341, -12.8774}, {-76.4341, -2.87735}, {-72.4341, -6.87735}, {-62.4341, 13.1226}, {-50.4341, -26.8774}, {-46.4341, -20.8774}, {-38.4341, -26.8774}, {-34.4341, -18.8774}, {-34.4341, 3.12265}, {-26.4341, 1.12265}, {-20.4341, 7.12265}, {-12.4341, 9.12265}, {-8.43408, 19.1226}, {1.56592, -4.87735}, {7.56592, -24.8774}, {19.5659, -6.87735}, {21.5659, 9.12265}, {31.5659, 13.1226}, {39.5659, -0.87735}, {43.5659, 11.1226}, {55.5659, 15.1226}, {63.5659, 27.1226}, {79.5659, -22.8774}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Rectangle(origin = {100, 0}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.05));
      end PiersonMoskowitzWave;

      model BrettschneiderWave
        /*  Brettschneider Spectrum - Wave elevation profile and excitation force */
        
        extends Modelica.Blocks.Icons.Block;
        import Modelica.Math.Vectors;
        OceanEngineeringToolbox.Internal.Connectors.WaveOutConn wconn;
        
        /* Environmental constants */
        constant Real pi = Modelica.Constants.pi "Mathematical constant pi";
        constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
        
        /*  Variable declarations */
        parameter String fileName;
        
        Real F_excRe[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.FexcRe", 1, 260));
        Real F_excIm[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.FexcIm", 1, 260));
        Real w[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.w", 1, 260));
        
        parameter Modelica.Units.SI.Length d = 100 "Water depth";
        parameter Modelica.Units.SI.Density rho = 1025 "Density of seawater";
        parameter Modelica.Units.SI.Length Hs = 2.5 "Significant Wave Height";
        parameter Modelica.Units.SI.AngularFrequency omega_min = 0.03141 "Lowest frequency component/frequency interval";
        parameter Modelica.Units.SI.AngularFrequency omega_max = 3.141 "Highest frequency component";
        parameter Modelica.Units.SI.AngularFrequency omega_peak = 0.9423 "Peak spectral frequency";
        parameter Integer n_omega = 100 "Number of frequency components";
        parameter Integer localSeed = 614657 "Local random seed";
        parameter Integer globalSeed = 30020 "Global random seed";
        parameter Real rnd_shft[n_omega] = OceanEngineeringToolbox.Internal.Functions.randomNumberGen(localSeed, globalSeed, n_omega);
        parameter Integer localSeed1 = 614757 "Local rand seed";
        parameter Integer globalSeed1 = 40020 "Global rand seed";
        parameter Real epsilon[n_omega] = OceanEngineeringToolbox.Internal.Functions.randomNumberGen(localSeed1, globalSeed1, n_omega) "Wave components phase shift";
        parameter Real Trmp = 200 "Interval for ramping up of waves during start phase";
        parameter Real omega[n_omega] = OceanEngineeringToolbox.Internal.Functions.frequencySelector(omega_min, omega_max, rnd_shft);
        parameter Real S[n_omega] = OceanEngineeringToolbox.Internal.Functions.spectrumGenerator_BRT(Hs, omega, omega_peak) "Spectral values for frequency components";
        parameter Modelica.Units.SI.Length zeta[n_omega] = sqrt(2*S*omega_min) "Wave amplitude component";
        parameter Real Tp[n_omega] = 2*pi./omega "Wave period components";
        parameter Real k[n_omega] = OceanEngineeringToolbox.Internal.Functions.waveNumber(d, omega) "Wave number component";
        Real ExcCoeffRe[n_omega] "Real component of excitation coefficient for frequency components";
        Real ExcCoeffIm[n_omega] "Imaginary component of excitation coefficient for frequency components";
        Real zeta_rmp[n_omega] "Ramp value for surface elevation";
        Modelica.Units.SI.Length SSE "Sea surface elevation";
        Modelica.Units.SI.Length SSE_unramp "Unramped sea surface elevation";
        
      equation
        /* Define amplitude for each frequency component - ramp time */
        for i in 1:n_omega loop
          if time < Trmp then
            zeta_rmp[i] = sin(pi/2*time/Trmp)*zeta[i];
          else
            zeta_rmp[i] = zeta[i];
          end if;

          /* Interpolate excitation coefficients (Re & Im) for each frequency component */
          ExcCoeffRe[i] = Modelica.Math.Vectors.interpolate(w2, F_excRe, omega[i]);
          ExcCoeffIm[i] = Modelica.Math.Vectors.interpolate(w2, F_excIm, omega[i]);
        end for;
        
        /* Define wave elevation profile (SSE) and excitation force */
        SSE = sum(zeta_rmp.*cos(omega*time - 2*pi*epsilon));
        SSE_unramp = sum(zeta.*cos(omega*time - 2*pi*epsilon));
        wconn.F_exc = sum((ExcCoeffRe.*zeta_rmp.*cos(omega*time - 2*pi*epsilon)) - (ExcCoeffIm.*zeta_rmp.*sin(omega*time - 2*pi*epsilon)))*rho*g;
        
        annotation(
          Icon(graphics = {Line(origin = {-50.91, 48.08}, points = {{-33.2809, -22.5599}, {-21.2809, -20.5599}, {-13.2809, 27.4401}, {6.71907, -20.5599}, {24.7191, -24.5599}, {42.7191, -24.5599}, {44.7191, -24.5599}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(origin = {-37, 51}, points = {{-51, 29}, {-51, -29}, {37, -29}}), Text(origin = {6, 55}, extent = {{-40, 17}, {40, -17}}, textString = "Hs"), Line(origin = {22, 4}, points = {{0, 22}, {0, -22}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-7.57, -61.12}, points = {{-82.4341, -12.8774}, {-76.4341, -2.87735}, {-72.4341, -6.87735}, {-62.4341, 13.1226}, {-50.4341, -26.8774}, {-46.4341, -20.8774}, {-38.4341, -26.8774}, {-34.4341, -18.8774}, {-34.4341, 3.12265}, {-26.4341, 1.12265}, {-20.4341, 7.12265}, {-12.4341, 9.12265}, {-8.43408, 19.1226}, {1.56592, -4.87735}, {7.56592, -24.8774}, {19.5659, -6.87735}, {21.5659, 9.12265}, {31.5659, 13.1226}, {39.5659, -0.87735}, {43.5659, 11.1226}, {55.5659, 15.1226}, {63.5659, 27.1226}, {79.5659, -22.8774}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Rectangle(origin = {100, 0}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.05));
      end BrettschneiderWave;

      model JonswapWave
        /*  JONSWAP Spectrum - Wave elevation profile and excitation force */
        
        extends Modelica.Blocks.Icons.Block;
        import Modelica.Math.Vectors;
        OceanEngineeringToolbox.Internal.Connectors.WaveOutConn wconn;
        
        /* Environmental constants */
        constant Real pi = Modelica.Constants.pi "Mathematical constant pi";
        constant Real g = Modelica.Constants.g_n "Acceleration due to gravity";
        
        /*  Variable declarations */
        parameter String fileName;
        
        Real F_excRe[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.FexcRe", 1, 260));
        Real F_excIm[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.FexcIm", 1, 260));
        Real w[:] = vector(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.w", 1, 260));
        
        /* Parameters */
        parameter Modelica.Units.SI.Length d = 100 "Water depth";
        parameter Modelica.Units.SI.Density rho = 1025 "Density of seawater";
        parameter Modelica.Units.SI.Length Hs = 2.5 "Significant Wave Height";
        parameter Modelica.Units.SI.AngularFrequency omega_min = 0.03141 "Lowest frequency component/frequency interval";
        parameter Modelica.Units.SI.AngularFrequency omega_max = 3.141 "Highest frequency component";
        parameter Modelica.Units.SI.AngularFrequency omega_peak = 0.9423 "Peak spectral frequency";
        parameter Real spectralWidth_min = 0.07 "Lower spectral bound for JONSWAP";
        parameter Real spectralWidth_max = 0.09 "Upper spectral bound for JONSWAP";
        parameter Integer n_omega = 100 "Number of frequency components";
        parameter Integer localSeed = 614657 "Local random seed";
        parameter Integer globalSeed = 30020 "Global random seed";
        parameter Real rnd_shft[n_omega] = OceanEngineeringToolbox.Internal.Functions.randomNumberGen(localSeed, globalSeed, n_omega);
        parameter Integer localSeed1 = 614757 "Local rand seed";
        parameter Integer globalSeed1 = 40020 "Global rand seed";
        parameter Real epsilon[n_omega] = OceanEngineeringToolbox.Internal.Functions.randomNumberGen(localSeed1, globalSeed1, n_omega) "Wave components phase shift";
        parameter Real Trmp = 200 "Interval for ramping up of waves during start phase";
        parameter Real omega[n_omega] = OceanEngineeringToolbox.Internal.Functions.frequencySelector(omega_min, omega_max, rnd_shft);
        parameter Real S[n_omega] = OceanEngineeringToolbox.Internal.Functions.spectrumGenerator_JONSWAP(Hs, omega, omega_peak, spectralWidth_min, spectralWidth_max) "Spectral values for frequency components";
        parameter Modelica.Units.SI.Length zeta[n_omega] = sqrt(2*S*omega_min) "Wave amplitude component";
        parameter Real Tp[n_omega] = 2*pi./omega "Wave period components";
        parameter Real k[n_omega] = OceanEngineeringToolbox.Internal.Functions.waveNumber(d, omega) "Wave number component";
        Real ExcCoeffRe[n_omega] "Real component of excitation coefficient for frequency components";
        Real ExcCoeffIm[n_omega] "Imaginary component of excitation coefficient for frequency components";
        Real zeta_rmp[n_omega] "Ramp value for surface elevation";
        Modelica.Units.SI.Length SSE "Sea surface elevation";
        Modelica.Units.SI.Length SSE_unramp "Unramped sea surface elevation";
        
      equation

        /* Define amplitude for each frequency component - ramp time */
        for i in 1:n_omega loop
          if time < Trmp then
            zeta_rmp[i] = sin(pi/2*time/Trmp)*zeta[i];
          else
            zeta_rmp[i] = zeta[i];
          end if;
          
          /* Interpolate excitation coefficients (Re & Im) for each frequency component */
          ExcCoeffRe[i] = Modelica.Math.Vectors.interpolate(w2, F_excRe, omega[i]);
          ExcCoeffIm[i] = Modelica.Math.Vectors.interpolate(w2, F_excIm, omega[i]);
        end for;
        
        /* Define wave elevation profile (SSE) and excitation force */
        SSE = sum(zeta_rmp.*cos(omega*time - 2*pi*epsilon));
        SSE_unramp = sum(zeta.*cos(omega*time - 2*pi*epsilon));
        wconn.F_exc = sum((ExcCoeffRe.*zeta_rmp.*cos(omega*time - 2*pi*epsilon)) - (ExcCoeffIm.*zeta_rmp.*sin(omega*time - 2*pi*epsilon)))*rho*g;
        
        annotation(
          Icon(graphics = {Line(origin = {-50.91, 48.08}, points = {{-33.2809, -22.5599}, {-21.2809, -20.5599}, {-13.2809, 27.4401}, {6.71907, -20.5599}, {24.7191, -24.5599}, {42.7191, -24.5599}, {44.7191, -24.5599}}, color = {255, 0, 0}, smooth = Smooth.Bezier), Line(origin = {-37, 51}, points = {{-51, 29}, {-51, -29}, {37, -29}}), Text(origin = {6, 55}, extent = {{-40, 17}, {40, -17}}, textString = "Hs"), Line(origin = {22, 4}, points = {{0, 22}, {0, -22}}, thickness = 1, arrow = {Arrow.None, Arrow.Filled}), Line(origin = {-7.57, -61.12}, points = {{-82.4341, -12.8774}, {-76.4341, -2.87735}, {-72.4341, -6.87735}, {-62.4341, 13.1226}, {-50.4341, -26.8774}, {-46.4341, -20.8774}, {-38.4341, -26.8774}, {-34.4341, -18.8774}, {-34.4341, 3.12265}, {-26.4341, 1.12265}, {-20.4341, 7.12265}, {-12.4341, 9.12265}, {-8.43408, 19.1226}, {1.56592, -4.87735}, {7.56592, -24.8774}, {19.5659, -6.87735}, {21.5659, 9.12265}, {31.5659, 13.1226}, {39.5659, -0.87735}, {43.5659, 11.1226}, {55.5659, 15.1226}, {63.5659, 27.1226}, {79.5659, -22.8774}}, color = {0, 0, 255}, smooth = Smooth.Bezier), Rectangle(origin = {100, 0}, fillColor = {85, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-20, 20}, {20, -20}})}, coordinateSystem(initialScale = 0.1)),
          experiment(StartTime = 0, StopTime = 400, Tolerance = 1e-06, Interval = 0.05));
      end JonswapWave;
      
    end IrregularWave;
  
  end WaveProfile;

  package Structures

    model WEC
      /* Solve Cummins' equation using state-space model of the radiation convolution integral */
      
      extends Modelica.Blocks.Icons.Block;
      OceanEngineeringToolbox.Internal.Connectors.WaveInConn wconn "Connector for wave elevation and excitation force" annotation(
        Placement(transformation(extent = {{-90, -10}, {-110, 10}})));
      OceanEngineeringToolbox.Internal.Connectors.DataCollector conn "Connector for velocity and radiation force" annotation(
        Placement(transformation(extent = {{-10, 90}, {10, 110}})));
      
      /* Parameters & variables */
      parameter String fileName;
      parameter Modelica.Units.SI.Mass M = scalar(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.m33", 1, 1)) "Total mass of the body (including ballast)";
      parameter Modelica.Units.SI.Mass Ainf = scalar(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.Ainf33", 1, 1)) "Added mass at maximum (cut-off) frequency";
      parameter Modelica.Units.SI.TranslationalSpringConstant Khs = scalar(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.Khs33", 1, 1)) "Hydrostatic stiffness";
      
      parameter Real A1[2, 2] = Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.ss_rad33.A", 2, 2) "State matrix";
      parameter Real B1[2, 1] = Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.ss_rad33.B", 2, 1) "Input matrix";
      parameter Real C1[1, 2] = Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.ss_rad33.C", 1, 2) "Output matrix";
      parameter Real D1 = scalar(Modelica.Utilities.Streams.readRealMatrix(fileName, "hydroCoeff.ss_rad33.D", 1, 1)) "Feed matrix";
      Real x[2,1] "State-space intermediate variables";
      Modelica.Units.SI.Length z "Heave displacement";
      Modelica.Units.SI.Velocity v_z "Heave velocity";
      Modelica.Units.SI.Acceleration a_z "Heave acceleration";
      Modelica.Units.SI.Force F_rad "Radiation Force";
      
    initial equation
      /* Define body at rest initially */
      z = 0;
      v_z = 0;
      
    equation
      v_z = der(z) "Heave velocity";
      a_z = der(v_z) "Heave acceleration";
      
      /* Radiation force state-space model */
      der(x) = (A1*x) + (B1*v_z);
      F_rad = scalar(C1*x) + (D1*v_z);
  
      /* Assemble Cummins' equation */
      ((M + Ainf)*a_z) + F_rad + (Khs*z) = wconn.F_exc;
  
      /* Connector declarations */
      conn.F_rad = F_rad;
      conn.v_z = v_z;
    end WEC;

  end Structures;

  package Internal
  
    package Functions
      /* Package defining explicit library functions */
  
      function waveNumber
        /* Function to iteratively compute the wave number from the frequency components */
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
  
      function randomNumberGen
        /* Function to generate random numbers from local and global seeds using XOR shift */
        input Integer ls = 614657 "Local seed";
        input Integer gs = 30020 "Global seed";
        input Integer n = 100 "Number of frequency components";
        output Real r64[n] "Random number vector";
      protected
        Integer state64[2](each start = 0, each fixed = true);
      algorithm
        state64[1] := 0;
        state64[2] := 0;
        for i in 1:n loop
          if i == 1 then
            state64 := Modelica.Math.Random.Generators.Xorshift64star.initialState(ls, gs);
            r64[i] := 0;
          else
            (r64[i], state64) := Modelica.Math.Random.Generators.Xorshift64star.random((state64));
          end if;
        end for;
      end randomNumberGen;
  
      function frequencySelector
        /* Function to randomly select frequency components */
        input Real omega_min "Frequency minima";
        input Real omega_max "Frequency maxima";
        input Real epsilon[:] "Random phase vector";
        output Real omega[size(epsilon, 1)] "Output vector of frequency components";
      protected
        parameter Real ref_omega[size(epsilon, 1)] = omega_min:(omega_max - omega_min)/(size(epsilon, 1) - 1):omega_max;
      algorithm
        omega[1] := omega_min;
        for i in 2:size(epsilon, 1) - 1 loop
          omega[i] := ref_omega[i] + epsilon[i]*omega_min;
        end for;
        omega[size(epsilon, 1)] := omega_max;
      end frequencySelector;
  
      function spectrumGenerator_PM
        /* Function to generate Pierson Moskowitz spectrum */
        input Real Hs = 1 "Significant wave height";
        input Real omega[:] "Frequency components";
        output Real spec[size(omega, 1)] "Spectral values for input frequencies";
      protected
        constant Real pi = Modelica.Constants.pi;
        constant Real g = Modelica.Constants.g_n;
      algorithm
        for i in 1:size(omega, 1) loop
          spec[i] := 0.0081*g^2/omega[i]^5*exp(-0.0358*(g/(Hs*omega[i]^2))^2);
        end for;
      end spectrumGenerator_PM;
  
      function spectrumGenerator_BRT
        /* Function to generate Brettschneider spectrum */
        input Real Hs = 1 "Significant wave height";
        input Real omega[:] "Frequency components";
        input Real omega_peak = 0.9423 "Peak spectral frequency";
        output Real spec[size(omega, 1)] "Spectral values for input frequencies";
      protected
        constant Real pi = Modelica.Constants.pi;
        constant Real g = Modelica.Constants.g_n;
      algorithm
        for i in 1:size(omega, 1) loop
          spec[i] := 1.9635*Hs^2*omega_peak^4/omega[i]^5*exp(-1.25*((omega_peak/omega[i])^4));
        end for;
      end spectrumGenerator_BRT;
  
      function spectrumGenerator_JONSWAP
        /* Function to generate JONSWAP spectrum */
        input Real Hs = 1 "Significant wave height";
        input Real omega[:] "Frequency components";
        input Real omega_peak = 0.9423 "Peak spectral frequency";
        input Real spectralWidth_min "Spectral width lower bound";
        input Real spectralWidth_max "Spectral width upper bound";
        output Real spec[size(omega, 1)] "Spectral values for input frequencies";
      protected
        constant Real pi = Modelica.Constants.pi;
        constant Real g = Modelica.Constants.g_n;
        constant Real gamma = 3.3;
        Real sigma;
        Real b;
      algorithm
        for i in 1:size(omega, 1) loop
          if omega[i] > omega_peak then
            sigma := spectralWidth_max;
          else
            sigma := spectralWidth_min;
          end if;
          b := exp(-0.5*(((omega[i] - omega_peak)/(sigma*omega_peak))^2));
          spec[i] := 0.0081*g^2/omega[i]^5*exp(-1.25*((omega_peak/omega[i])^4))*gamma^b;
        end for;
      end spectrumGenerator_JONSWAP;
    end Functions;
  
    package Connectors
      /* Package defining library connectors between models */
  
      connector WaveOutConn
        /* Output datastream - wave elevation & excitation force */
        Modelica.Blocks.Interfaces.RealOutput F_exc;
      end WaveOutConn;
  
      connector WaveInConn
        /* Input datastream - wave elevation & excitation force */
        Modelica.Blocks.Interfaces.RealInput F_exc;
      end WaveInConn;
  
      connector DataCollector
        /* Output datastream - velocity and radiation force */
        Modelica.Blocks.Interfaces.RealOutput F_rad;
        Modelica.Blocks.Interfaces.RealOutput v_z;
      end DataCollector;
    end Connectors;

  end Internal;

  package Tutorial
    /* Sample simulation models */
    
    model sample1
      /* Single body, regular waves */
      
      parameter String filePath = "C:/.../hydroCoeff.mat";
      OceanEngineeringToolbox.WaveProfile.RegularWave.LinearWave Reg1(fileName = filePath, Hs = 2.5, Trmp = 50);
      OceanEngineeringToolbox.Structures.WEC WEC1(fileName = filePath);
    equation
      connect(Reg1.wconn.F_exc, WEC1.wconn.F_exc);
      annotation(
        Line(points = {{-40, 30}, {-20, 30}, {-20, 0}, {40, 0}, {40, 0}}),
        experiment(StartTime = 0, StopTime = 500, Tolerance = 1e-06, Interval = 0.1));
    end sample1;
    
    model sample2
      /* Single body, irregular waves with PM spectrum */
      
      parameter String filePath = "C:/.../hydroCoeff.mat";
      OceanEngineeringToolbox.WaveProfile.IrregularWave.PiersonMoskowitzWave PM1(fileName = filePath, Hs = 2.5, n_omega = 100, Trmp = 50);
      OceanEngineeringToolbox.Structures.WEC WEC1(fileName = filePath);
    equation
      connect(PM1.wconn.F_exc, WEC1.wconn.F_exc);
      annotation(
        Line(points = {{-40, 30}, {-20, 30}, {-20, 0}, {40, 0}, {40, 0}}),
        experiment(StartTime = 0, StopTime = 200, Tolerance = 1e-06, Interval = 0.1));
    end sample2;
    
  end Tutorial;

  package Simulations
    /* Directory for user-defined simulation models */
  end Simulations;
end OceanEngineeringToolbox;
/*  Modelica Ocean Engineering Toolbox (OET)
    Developed at:
          Sys-MoDEL, 
          University of New Brunswick, Fredericton
          New Brunswick, E3B 5A3, Canada
    Copyright under the terms of the GNU General Public License
*/
