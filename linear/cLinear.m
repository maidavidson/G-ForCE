function [c]=cLinear(T,d)
%% Function or program name: cLinear
% CLINEAR G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       T — input argument (description/units to be verified)
%       d — input argument (description/units to be verified)
%
%   Outputs:
%       c — output argument (description/units to be verified)
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

%CLINEAR G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% [c]=cLinear(T,d);
% Approximate solution to the linear wave dispersion equation after Hunt 1979
% Mark Davidson
g=9.81; % Gravity
k=kLinear(T,d); % Wave number
c=(g*T./(2*pi)).*tanh(k.*d);
