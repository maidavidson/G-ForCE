%% Function below...
function [y]=myInterp(xo,yo,x)
%% Function or program name: myInterp
% Function below...
%
%   Inputs:
%       xo — input argument (description/units to be verified)
%       yo — input argument (description/units to be verified)
%       x — input argument (description/units to be verified)
%
%   Outputs:
%       y — output argument (description/units to be verified)
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

%MYINTERP G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
if isempty(xo) || isempty (yo); y=[];return;end

    y=(yo(1)*(xo(2)-x)+yo(2)*(x-xo(1)))/(xo(2)-xo(1));
end
