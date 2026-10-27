%% DOC calculation
function [DOC,H12,T12]=Hallermeier(H,T,dt)
%% Function or program name: Hallermeier
% DOC calculation
%
%   Inputs:
%       H — input argument (description/units to be verified)
%       T — input argument (description/units to be verified)
%       dt — input argument (description/units to be verified)
%
%   Outputs:
%       DOC — output argument (description/units to be verified)
%       H12 — output argument (description/units to be verified)
%       T12 — output argument (description/units to be verified)
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

%HALLERMEIER G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% deep water H array [m]
% deep water T array [s]
% wave sampling interval (must be regular)
% DOC - depth of closure [m]
% H12 - 12 hour in one year wave height [m]
% T12 - same for T [s]

H=sort(H,'descend');
T=sort(T,'descend');
Npts=length(H);
Ndays=(Npts*dt);
indexVal=round(Ndays*0.5/365.25);
indexVal=max(indexVal,1);
H12=H(indexVal);
T12=T(indexVal);
DOC= 2.28 * H12 - 68.5 * (H12.^2/ (9.81*T12.^2));

end
