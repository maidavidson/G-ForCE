function [TC]=transportCost(x,y1,y2,ym)
%% Function or program name: TransportCost
% TRANSPORTCOST G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       x — input argument (description/units to be verified)
%       y1 — input argument (description/units to be verified)
%       y2 — input argument (description/units to be verified)
%       ym — input argument (description/units to be verified)
%
%   Outputs:
%       TC — output argument (description/units to be verified)
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

%TRANSPORTCOST G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.

% x  = chainage
% y1 = measured initial profile
% y2 = measured final profile
% ym = modelled final profile
% TC = dimensionless transport cost skill score

dx=x(2)-x(1);
boundaryType=   [1 1]; % 0 = closed (Neumann), 1 = free (Dirichlet)
boundaryValues= [0 0]; % NB Specify a gradient for Neumann and a constant for Dirichlet

% Transport cost between measured profile
[phi]=Poisson(dx,y1,y2,boundaryType,boundaryValues);
 q1=gradient(phi,x);

 % Transport cost between measured profile
 [phi]=Poisson(dx,ym,y2,boundaryType,boundaryValues);
 q2=gradient(phi,x);

 TC=1-(rms(q2)/rms(q1));

return
end
%%
%[text] ## Poisson Equation

function [phi]=Poisson(dx,y1,y2,boundaryType,boundaryValues)
%POISSON G_ForCE support routine; see function signature for inputs and outputs.
% dx - x resolution
% y1 - Model
% y2 - Measurement
% boundaryType: 0 - closed (Neumann), 1 - Free (Dirichlet), e.g. boundaryType=[0 1] - [closed
% free]
% boundaryValue - Specify gradient for Neumann or constant value for Dirichlet

% Source function
S=dx.^2*(y1-y2);

% Tridiagonal matrix coefficients a, b and c
a=ones(size(S)); b=-2*a; c=a;

% Define left and right hand boundary conditions
% l.h.s. boundary
if boundaryType(1)==0
    % Closed (Neumann) boundary conditions
    c(1)=2; S(1)=S(1)+2*dx*boundaryValues(1);
else
    % Free (Dirichlet) boundary conditions
    S(1)=S(1)-boundaryValues(1);
end

% r.h.s. boundary
if boundaryType(2)==0
    % Closed (Neumann) boundary conditions
    a(end)=2; S(end)=S(end)+2*dx*boundaryValues(2);
else
    % Free (Dirichlet) boundary conditions
    S(end)=S(end)-boundaryValues(2);
end

% Assemble tridiagonal matrix and derive implicit solution
[phi]=implicitMD2(a,b,c,S);

end % End subroutine

function [y]=implicitMD2(a,b,c,S)
%IMPLICITMD2 G_ForCE support routine; see function signature for inputs and outputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implicit solution to tridiagonal equations                            %
% Input:                                                                %                                                  %
% a, b, c implicit inversion coefficients.                              %
% S: Source function f(x).                                              %
% Output:                                                               %
% y are approximate solution at grid points                             %
% Where: a*phi(i-1)+b*phi(i)+c*phi(i+1)=S                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S=S(:);
n = length(S);
A = sparse(n,n);
y = zeros(1,n);

A(1,[1 2])       = [b(1) c(1)]; % lhs boundary
for i=2:n-1
    A(i,[i-1 i i+1]) = [a(i) b(i) c(i)];
end
A(n,[n-1 n])     = [a(n) b(n)]; % rhs boundary

y = A\S;
end
%[text]
%[text]

%[appendix]{"version":"1.0"}
%---
