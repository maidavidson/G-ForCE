function [Ir]=IrribarenNo(Ho,Tp,tanb)
%% Function or program name: IrribarenNo
% IRRIBARENNO G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       Ho — input argument (description/units to be verified)
%       Tp — input argument (description/units to be verified)
%       tanb — input argument (description/units to be verified)
%
%   Outputs:
%       Ir — output argument (description/units to be verified)
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

%IRRIBARENNO G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
if Ho==0;Ir=inf;return;end
    [Lo]=LoLinear(Tp/1.1);
    Ir=tanb/sqrt(Ho/Lo);
end
