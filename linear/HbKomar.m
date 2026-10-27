% [Hb]=HbKomar(Ho,T);
% Hb = breaking wave height (m)
% Ho = offshore significant wave height (m)
% T  = peak wave period (s)
% P.D. Komar (1974) equation
function [Hb]=HbKomar(Ho,T);
%% Function or program name: HbKomar
% [Hb]=HbKomar(Ho,T);
%
%   Inputs:
%       Ho — input argument (description/units to be verified)
%       T — input argument (description/units to be verified)
%
%   Outputs:
%       Hb — output argument (description/units to be verified)
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

%HBKOMAR G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
Hb = 0.39*9.81.^0.2*(T.*Ho.^2).^0.4;
