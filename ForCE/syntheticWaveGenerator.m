%% Function or program name: syntheticWaveGenerator
% Synthetic wave generator (Mark Davidson 2025)
%
%   Inputs:
%       None (script; see original code for workspace/file dependencies)
%
%   Outputs:
%       None (no formal return arguments)
%
%   Author:  Dr Mark Davidson
%   Project: NERC Gravel Beach Project
%   Contact: mdavidson@plymouth.ac.uk
%   Date:    26-Oct-2026
%   Version: 1.0.0
%
%   Toolbox Requirements:
%       Not independently verified; consult called functions and MATLAB dependencies.
%
%   See also: G_ForCE_v1_1, runG_ForCE
%
%==========================================================================
% REVISION HISTORY
%==========================================================================
% Version | Date        | Author       | Description
% --------|-------------|--------------|-----------------------------------
% 1.0.0   | 26-Oct-2026 | M Davidson   | Standardised documentation header.
%==========================================================================
% end header

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
%% Synthetic wave generator (Mark Davidson 2025)
disp('Synthetic wave generator')

%% Get waves file
 disp('Please load the example wave file containg the structure wave with fields H, Tp and Dir')
 [file,path] = uigetfile('*.mat');

% Setup paths. Modifiy first one with GUI later
 addpath('./ForCE/')

% Load data
 load(fullfile(path,file))

% Define the number of years in the measured/modelled wave data
endDate=wave.dates(end);
poolSize=ceil((endDate-wave.dates(1))/365);

% Call the synthetic wave making routine
[forecastDates, Var1Forecast,Var2Forecast,Var3Forecast]=makingWavesFunction(wave.dates,wave.H,wave.Tp,wave.Dir,wave.dates(end),poolSize,100*365,wave.dt);

% refresh wave structure to include forecast values.
wave.dates=forecastDates;
wave.H=Var1Forecast;
wave.Tp=Var2Forecast;
wave.Dir=Var3Forecast;
% flag synthetic data
wave.flag=zeros(size(wave.dates));
jj=find(wave.dates>endDate);
wave.flag(jj)=1;

disp('Synthetic waves file generated with 100 years of synthetic wave tagged onto the original date file')
disp('now save the file in the input directory for your site')
uisave('wave')
