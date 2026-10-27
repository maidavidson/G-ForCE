%% Find solution via matrix inversion
function [y]=implicit(a,b,c,S)
%% Function or program name: implicit
% Find solution via matrix inversion
%
%   Inputs:
%       a — input argument (description/units to be verified)
%       b — input argument (description/units to be verified)
%       c — input argument (description/units to be verified)
%       S — input argument (description/units to be verified)
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


% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implicit solution to tridiagonal equations
% Input:
% a, b, c implicit inversion coefficients.
% S: Source function f(x).
% Output:
% y are approximate solution at grid points
% Where: a*y(i-1)+b*y(i)+c*y(i+1)=S
% Written by Mark Davidson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [n,m] = size(S);
 A = sparse(n,n);

 A(1,[1 2])       = [b(1) c(1)]; % lhs boundary
for i=2:n-1
 A(i,[i-1 i i+1]) = [a(i) b(i) c(i)];
end
 A(n,[n-1 n])     = [a(n) b(n)]; % rhs boundary

 y = A\S;
end
