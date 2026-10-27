function [BSS]=brierss(xm,x,xb,dx)
%% Function or program name: brierss
% Brier Skils Function
%
%   Inputs:
%       xm — input argument (description/units to be verified)
%       x — input argument (description/units to be verified)
%       xb — input argument (description/units to be verified)
%       dx — input argument (description/units to be verified)
%
%   Outputs:
%       BSS — output argument (description/units to be verified)
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
% Brier Skils Function

% xm - measured data
% x - model data
% xb - benchmrk data
% dx - measurement error

BSS = 1 - ...
    ( ...
    mean( ( ( abs(x-xm) - dx ) ).^2) /...
    mean( (xb-xm).^2) ...
    );
