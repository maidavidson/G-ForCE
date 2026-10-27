%% Avelanching algorithm
% Author Mark Davidson
function [zout,iterations]=avelanche(x,zin,model)
%% Function or program name: avelanche
% Avelanching algorithm
%
%   Inputs:
%       x — input argument (description/units to be verified)
%       zin — input argument (description/units to be verified)
%       model — input argument (description/units to be verified)
%
%   Outputs:
%       zout — output argument (description/units to be verified)
%       iterations — output argument (description/units to be verified)
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

%AVELANCHE G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.

% Initialise variables
 kDiff=model.dx.^2; % Diffussion coefficient
 dx=model.dx;
 maxIteration=model.maxIteration; % Max no. of iterations
 distSmoothed=model.distSmoothed; % Region either side of gradient excess to be considerred
 dt=1;
 fact=0.5;
 zout=zin; % Smoothed profile array
 ii=1; % interation counter

%% Loop twise here 1) For dune and 2) Active beach
  maxGrad=model.maxGrad;
  iActive=1:length(zin);

% Compute profile gradient and find where max gradient is exceeded
 profGrad=gradient(zin,x);
 jj=find(abs(profGrad(iActive))>maxGrad);

% Main iterative loop (Iterate until smoothed to model.maxGrad)
 while  length(jj)>1 && ii<maxIteration

% Warn user of avelanche
  disp('************* Avelanche **************')

% Define region around excess gradient
  npts=ceil(distSmoothed/model.dx);
  istart=1;
  ifin=max(jj)+npts;
  istart=max(istart,2);ifin=min(ifin,length(zin)-5);

% Define points to be avelached
  hPrime=zout(istart:ifin);

% Break if no data
 if length(hPrime)<5;break;end

% Crank-Nicholson Solution Diffusion
% Define boundary type and values:
% 0 - closed
% 1 - free
%  boundaryType=[1 1]; % 0 = closed (Neumann), 1 = free (Dirichlet)
%  boundaryValues=[zout(istart-1) zout(ifin+1)]; % NB Specify a gradient for Neumann and a constant for Dirichlet
 boundaryType=[0 0]; % 0 = closed (Neumann), 1 = free (Dirichlet)
 boundaryValues=[profGrad(istart) profGrad(ifin)]; % NB Specify a gradient for Neumann and a constant for Dirichlet
 [hPrime]=crankNicholson(hPrime,kDiff,dt,dx,fact,boundaryType,boundaryValues);

% Insert smoothed profile into original
  zout(istart:ifin)=hPrime;

% recompute new gradient and region in excess of maxGrad
  clear jj profGrad
  profGrad=gradient(zout,x);
   jj=find(abs(profGrad(iActive))>maxGrad);

% Iteration counter
  ii=ii+1;

 end % end while loop
%% iteration counter
 iterations=ii-1;

end
%%
