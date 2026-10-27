function [mean_bearing, delta_deg]=convertAngle(b,mean_bearing)
%% Function or program name: convertAngle
% CONVERTANGLE G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       b — input argument (description/units to be verified)
%       mean_bearing — input argument (description/units to be verified)
%
%   Outputs:
%       mean_bearing — output argument (description/units to be verified)
%       delta_deg — output argument (description/units to be verified)
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

%CONVERTANGLE G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.

if nargin==1

% convert compass bearings to mathematical angles (0 = +x, CCW positive)
theta = deg2rad(90 - b);      % math convention radians

% circular mean (complex)
z = mean(exp(1i*theta), 'omitnan');

% mean bearing back to compass convention (0..360)
mean_bearing = mod(90 - rad2deg(angle(z)), 360);

% signed shortest difference bearing - mean_bearing in (-180,180]
% compute as circular difference in math radians then convert to degrees
diff_rad = angle(exp(1i*(theta - angle(z))));    % gives value in (-pi,pi]
delta_deg = - rad2deg(diff_rad);                   % corresponds to b - mean_bearing


else

d=b-mean_bearing;
delta_deg=mod(d+180,360) - 180;

end
