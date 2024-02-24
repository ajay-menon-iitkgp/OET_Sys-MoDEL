%% Simulation Data
simu = simulationClass();                       % Initialize simulationClass
simu.simMechanicsFile = 'floatLineAttenuator.slx';            % Simulink Model File
simu.startTime = 0;                             % Simulation Start Time [s]
simu.rampTime = 100;                            % Wave Ramp Time [s]
simu.endTime = 400;                             % Simulation End Time [s]   
simu.dt = 0.1;                                  % Simulation Time-Step [s]

%% Wave Information
% % Regular Waves 
% waves = waveClass('regular');                   % Initialize waveClass                                 
% waves.height = 2.5;                                  % Wave Height [m]
% waves.period = 8;                                    % Wave Period [s]

% Waves with imported wave elevation time-history  
waves = waveClass('elevationImport');          % Create the Wave Variable and Specify Type
waves.elevationFile = 'elevationData.mat';     % Name of User-Defined Time-Series File [:,2] = [time, eta]

%% Body Data
% Float Body
body(1) = bodyClass('hydroData/rm3.h5');      % Initialize bodyClass for Flap 
body(1).geometryFile = 'geometry/float.stl';     % Geometry File
body(1).mass = 'equilibrium';                   
body(1).inertia = [20907301 21306090.66 37085481.11];  % Moment of Inertia [kg*m^2]

%% PTO and Constraint Parameters
% Floating (3DOF) Joint
constraint(1)= constraintClass('Constraint1');  % Initialize constraintClass for Constraint1
constraint(1).location = [0 0 0];                  % Constraint Location [m]