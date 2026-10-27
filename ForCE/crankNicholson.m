%% Crank-Nicholson finite-difference diffusion routine.
% [fux conservative + variable diffusion coefficient]. Written by Mark Davidson
 function [hPrime]=crankNicholson(hPrime,kDiff,dt,dx,fact,boundaryType,boundaryValues)
%% Function or program name: crankNicholson
% Crank-Nicholson finite-difference diffusion routine.
%
%   Inputs:
%       hPrime — input argument (description/units to be verified)
%       kDiff — input argument (description/units to be verified)
%       dt — input argument (description/units to be verified)
%       dx — input argument (description/units to be verified)
%       fact — input argument (description/units to be verified)
%       boundaryType — input argument (description/units to be verified)
%       boundaryValues — input argument (description/units to be verified)
%
%   Outputs:
%       hPrime — output argument (description/units to be verified)
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

%CRANKNICHOLSON G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.

% Find number of points and initialise temporary arrays (hPrimeNew)
 Nx=length(hPrime);
 hPrimeNew=zeros(size(hPrime));
 a=hPrime;b=a;c=a;
 a2=a;b2=a;c2=a;

% Define combined diffusion coefficient (Dcoeff)
  if length(kDiff)==1
   Dcoeff=(kDiff*dt./dx.^2).*ones(size(hPrime));
  else
   Dcoeff=kDiff*dt./dx.^2;
  end

% Define left and right hand boundary conditions for EXPLICIT components
% only
% l.h.s. boundary
 if boundaryType(1)==0
  % Closed (Neumann) boundary conditions
  a(1)=-2*dt/dx;
  b(1)=(1 - (1-fact)*Dcoeff(1) - (1-fact)*Dcoeff(2));
  c(1)=  2 * (1-fact)*Dcoeff(2);
else
 % Free (Dirichlet) boundary conditions
  a(1)=(1-fact)*Dcoeff(1);
  b(1)=(1 - (1-fact)*Dcoeff(1) - (1-fact)*Dcoeff(2));
  c(1)=(1-fact)*Dcoeff(2);
end
% explicit computation (start boundary)
 hPrimeNew(1)= ...
   a(1) * boundaryValues(1) +...
   b(1) * hPrime(1)                 +...
   c(1) * hPrime(2);

% r.h.s. boundary
 if boundaryType(2)==0
 % Closed (Neumann) boundary conditions
  a(end)=2*(1-fact)*Dcoeff(end-1);
  b(end)=(1 - (1-fact)*Dcoeff(end-1) - (1-fact)*Dcoeff(end));
  c(end)=2*dt/dx;
 else
 % Free (Dirichlet) boundary conditions
  a(end)=(1-fact)*Dcoeff(end-1);
  b(end)=(1 - (1-fact)*Dcoeff(end-1) - (1-fact)*Dcoeff(end));
  c(end)=(1-fact)*Dcoeff(end);
end
% explicit computation (end boundary)
 hPrimeNew(end)=...
   a(end) * hPrime(end-1) +...
   b(end) * hPrime(end)     +...
   c(end) * boundaryValues(2);

 % Compute non-boundary explicit terms and assemple implicit coefficients a2, b2 and c2
 % NB Could replace this with more efficient matrix multiplication omitting
 % loop
for jj=2:Nx-1; % Start EXPLICIT diffusion loop
 % Compute EXPLICIT diffusion coefficients
  a(jj)=(1-fact)*Dcoeff(jj);
  b(jj)=(1 - (1-fact)*Dcoeff(jj) - (1-fact)*Dcoeff(jj+1));
  c(jj)=(1-fact)*Dcoeff(jj+1);

 % Compute IMPLICT diffusion coefficients
  a2(jj)=-fact*Dcoeff(jj);
  b2(jj)=(1 + fact*Dcoeff(jj) + fact*Dcoeff(jj+1));
  c2(jj)=-fact*Dcoeff(jj+1);

 % EXPLICIT diffusion component (hPrimeNew)
  hPrimeNew(jj)=...
    a(jj) * hPrime(jj-1) +...
    b(jj) * hPrime(jj) +...
    c(jj) * hPrime(jj+1);
end % end EXPLICIT diffusion loop

% Define left and right hand boundary conditions for IMPLICIT components
% only
% l.h.s. boundary
 if boundaryType(1)==0
  % Closed (Neumann) boundary conditions
  a2(1)=0; % unused
  b2(1)=(1 + fact*Dcoeff(1) + fact*Dcoeff(2));
  c2(1)=-2*fact*Dcoeff(2);
  hPrimeNew(1)=hPrimeNew(1) - (2*fact*Dcoeff(1)*dx*boundaryValues(1));
else
 % Free (Dirichlet) boundary conditions
  a2(1)=0; % unused
  b2(1)=(1 + fact*Dcoeff(1) + fact*Dcoeff(2));
  c2(1)=-fact*Dcoeff(2);
  hPrimeNew(1)=hPrimeNew(1) + fact*Dcoeff(1)*boundaryValues(1);
 end

% r.h.s. boundary
 if boundaryType(2)==0
 % Closed (Neumann) boundary conditions
  a2(end)=-2*fact*Dcoeff(end-1);
  b2(end)=(1 + fact*Dcoeff(end-1) + fact*Dcoeff(end));
  c2(end)=0; % unused
  hPrimeNew(end)=hPrimeNew(end) + (fact*Dcoeff(end)*2*dx*boundaryValues(2));
 else
 % Free (Dirichlet) boundary conditions
  a2(end)=-fact*Dcoeff(end-1);
  b2(end)=(1 + fact*Dcoeff(end-1) + fact*Dcoeff(end));
  c2(end)=0; % unused
  hPrimeNew(end)=hPrimeNew(end) + fact*Dcoeff(end)*boundaryValues(2);
 end

 % SOLUTION - Invert tridiagonal matrix and update and return hPrime
  [hPrime] = implicit(a2,b2,c2,hPrimeNew');hPrime=hPrime';

 end % end function
