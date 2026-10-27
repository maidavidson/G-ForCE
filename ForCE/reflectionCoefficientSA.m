function [R]=reflectionCoefficientSA(Ir,a,b)
%% Function or program name: reflectionCoefficientSA
% REFLECTIONCOEFFICIENTSA G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       Ir — input argument (description/units to be verified)
%       a — input argument (description/units to be verified)
%       b — input argument (description/units to be verified)
%
%   Outputs:
%       R — output argument (description/units to be verified)
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

%REFLECTIONCOEFFICIENTSA G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
if isinf(Ir);R=a;return;end
Ir2=Ir.^2;
R=(a*Ir2)/(b+Ir2);

end
