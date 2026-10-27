function [P]=PLinear(E,C,n)
%% Function or program name: PLinear
% PLINEAR G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       E — input argument (description/units to be verified)
%       C — input argument (description/units to be verified)
%       n — input argument (description/units to be verified)
%
%   Outputs:
%       P — output argument (description/units to be verified)
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

%PLINEAR G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% [P]=PLinear(E,C,n);
% Linear Wave Power Function
% E - Wave energy, C - wave celerity, n - linear n-value
% Mark Davidson
P=E.*C.*n;
