function [k]=kLinear(T,d)
%% Function or program name: kLinear
% KLINEAR G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       T — input argument (description/units to be verified)
%       d — input argument (description/units to be verified)
%
%   Outputs:
%       k — output argument (description/units to be verified)
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

%KLINEAR G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% [k]=kLinear(T,d);
% Approximate solution to the linear wave dispersion equation after Hunt 1979
% Mark Davidson
g=9.81;     % gravity
y=4.03*d./(T.^2);
kd2=y.^2 + y ./ ...
    (1+(0.666.*y)+(0.355.*y.^2)+(0.161.*y.^3)+(0.0632.*y.^4)+(0.0218.*y.^5)+(0.00564.*y.^6));
k=sqrt(kd2)./d;
