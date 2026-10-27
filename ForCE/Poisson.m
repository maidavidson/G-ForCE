%% Solve Poisson Eq.
function [phi]=Poisson(dx,y1,y2,boundaryType,boundaryValues)
%% Function or program name: Poisson
% Solve Poisson Eq.
%
%   Inputs:
%       dx — input argument (description/units to be verified)
%       y1 — input argument (description/units to be verified)
%       y2 — input argument (description/units to be verified)
%       boundaryType — input argument (description/units to be verified)
%       boundaryValues — input argument (description/units to be verified)
%
%   Outputs:
%       phi — output argument (description/units to be verified)
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

%POISSON G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% dx - x resolution
% y1 - Model
% y2 - Measurement
% boundaryType: 0 - closed (Neumann), 1 - Free (Dirichlet), e.g. boundaryType=[0 1] - [closed
% free]
% boundaryValue - Specify gradient for Neumann or constant value for Dirichlet

% Source function
 S=dx.^2*(y1-y2);

 % Length of array
  [n,m]=size(S);

% Tridiagonal matrix coefficients a, b and c
 a=ones(1,n); b=-2*a; c=a;

% Define left and righ hand boundary conditions
% l.h.s. boundary
if boundaryType(1)==0
 % Closed (Neumann) boundary conditions
 c(1)=2; S(:,1)=S(:,1)+2*dx*boundaryValues(1);
else
 % Free (Dirichlet) boundary conditions
 S(:,1)=S(:,1)-boundaryValues(1);
end

% r.h.s. boundary
if boundaryType(2)==0
 % Closed (Neumann) boundary conditions
 a(end)=2; S(:,end)=S(:,end)+2*dx*boundaryValues(2);
else
 % Free (Dirichlet) boundary conditions
 S(:,end)=S(:,end)-boundaryValues(2);
end

% Assemble tridiagonal matrix and derive implicit solution
[phi]=implicit(a,b,c,S);

end % End subroutine
