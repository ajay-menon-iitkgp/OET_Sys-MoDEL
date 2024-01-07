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
