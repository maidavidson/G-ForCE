function [wi]=infiltrationVelocity(D50,C,hg,poreSpace)
%% Function or program name: infiltrationVelocity
% INFILTRATIONVELOCITY G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       D50 — input argument (description/units to be verified)
%       C — input argument (description/units to be verified)
%       hg — input argument (description/units to be verified)
%       poreSpace — input argument (description/units to be verified)
%
%   Outputs:
%       wi — output argument (description/units to be verified)
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

%INFILTRATIONVELOCITY G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% C=0.01-0.05 s^-1
% hydraulic gradient (hg) = 0.5-3
% poreSpace = 0.3-0.4 gravel
if isempty(C)
    C=0.025;
end

if isempty(hg)
    if D50<0.002
        hg=1;
    else
        hg=3;
    end
end

if isempty(poreSpace)
   if D50<0.002
    poreSpace=0.4;
   else
    poreSpace=0.35;
   end
end

    K=C*D50.^2;
    q=K*hg;
    wi=q/poreSpace;

    % Sanity check


end
