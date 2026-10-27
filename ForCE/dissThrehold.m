function [Dcrit]=dissThrehold(Tp, model,layerNo)
%% Function or program name: dissThrehold
% DISSTHREHOLD G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       Tp — input argument (description/units to be verified)
%       model — input argument (description/units to be verified)
%       layerNo — input argument (description/units to be verified)
%
%   Outputs:
%       Dcrit — output argument (description/units to be verified)
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

%DISSTHREHOLD G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% f - friction factor >=1
% D50 - median grain size [m]
% Tp - Peak wave period [s]
% Dcrit - Critical Dissipation for sediment transport

if nargin==2
    layerNo=1;
end

Dcrit= (2*model.f*(model.sedimentDensity-model.waterDensity)*9.81*model.D50(layerNo).^2)./...
        (3*Tp);

end
