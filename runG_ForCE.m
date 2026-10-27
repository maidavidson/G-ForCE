%% Function or program name: runG_ForCE
% Module to open a setup file and run G-ForCE
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
%       Not independently verified; consult called functions and MATLAB 
%       dependencies.
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

% This code was written by Mark Davidson and amended for application to 
% gravel and mixed beaches as part of the NERC Gravel Beach Project.
% Module to open a setup file and run G-ForCE
clear all
close all
% Get setup file
 disp('Please load the model setup file found in the project directory')
 [file,path] = uigetfile('*.m');
 if isequal(file,0) || isequal(path,0)
     disp('Setup file selection cancelled.');
     return
 end
% caffeinate
 system('caffeinate &');
% Evauate the setup file
 run(fullfile(path,file))
% Run G-ForCE
 run(filename.executable)
