function [runup]=runupSimple(Ho,c)
%% Function or program name: runupSimple
% RUNUPSIMPLE G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       Ho — input argument (description/units to be verified)
%       c — input argument (description/units to be verified)
%
%   Outputs:
%       runup — output argument (description/units to be verified)
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

%RUNUPSIMPLE G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
if nargin==1
    c=0.5;
end
    runup=c*Ho;
end
