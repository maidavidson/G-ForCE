function [Hd]=HoLinear(H1,T,d1,d);
%% Function or program name: HoLinear
% HOLINEAR G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       H1 — input argument (description/units to be verified)
%       T — input argument (description/units to be verified)
%       d1 — input argument (description/units to be verified)
%       d — input argument (description/units to be verified)
%
%   Outputs:
%       Hd — output argument (description/units to be verified)
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

%HOLINEAR G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% [Hd]=HoLinear(H,T,d1,d);
% Hd - wave height at depth d, H1 - Reference wave height at depth d1, T - Period, d - depth
% where the deep water wave height is required (enter inf for deep water)
% NB Assumes normal incidence over planar bed
% Mark Davidson
if d ==inf; d=10.^17;  end
if d1==inf; d1=10.^17; end
c1 =cLinear(T,d1);
cd =cLinear(T,d);
k1 =kLinear(T,d1);
kd =kLinear(T,d);
n1 =nLinear(k1,d1);
nd =nLinear(kd,d);
Hd=sqrt((n1.*c1) ./ (nd.*cd)) .* H1;
