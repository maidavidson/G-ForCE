function [w,wi]=cosexpdecay(N)
%% Function or program name: cosexpdecay
% COSEXPDECAY G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       N — input argument (description/units to be verified)
%
%   Outputs:
%       w — output argument (description/units to be verified)
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

%COSEXPDECAY G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% cos x exp Decay function
% error catch
if N<2;w=0;wi=1;return;end

% define wavelength
L=2 * (N - 1);
% wave number
k=2*pi/L;
% x values
x= 0:N-1;
% exp
y1=exp(-2*k*x);
%sin
y2=0.5*(cos(k*x)+1);
% weight
w=y1.*y2;
% inverse weight
wi=1-w;
end
