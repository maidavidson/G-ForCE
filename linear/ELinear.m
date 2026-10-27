function [E]=ELinear(H)
%% Function or program name: ELinear
% ELINEAR G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       H — input argument (description/units to be verified)
%
%   Outputs:
%       E — output argument (description/units to be verified)
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

%ELINEAR G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% [E]=ELinear(H);
% Linear Wave Energy Function
% H - wave height, E- Wave energy
% Mark Davidson
E=(1025*9.81*H.^2)./16;
