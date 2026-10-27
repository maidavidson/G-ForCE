% NOTE: This standalone implementation scales BOTH setup and swash by swashFact.
% The local StockdonHolman in WaveModelBJ_swash.m scales swash only.
% Do not consolidate without validating the intended physical formulation.
function [R2,setup,swash]=StockdonHolman(Ho,T,beta,swashFact)
%% Function or program name: StockdonHolman
% STOCKDONHOLMAN G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       Ho — input argument (description/units to be verified)
%       T — input argument (description/units to be verified)
%       beta — input argument (description/units to be verified)
%       swashFact — input argument (description/units to be verified)
%
%   Outputs:
%       R2 — output argument (description/units to be verified)
%       setup — output argument (description/units to be verified)
%       swash — output argument (description/units to be verified)
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

%STOCKDONHOLMAN G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
g=9.98;
Lo=(g*T.^2)/(2*pi);
setup=swashFact * 0.35 * beta * sqrt(Ho.*Lo);
swash=swashFact * 0.5 * sqrt(Ho .* Lo * (0.563 * beta.^2 + 0.004));
R2=(setup+swash);
end
