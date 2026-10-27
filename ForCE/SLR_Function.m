function [heq,Vold,model,dxSLR,SLR,SLROld,flagUpdate]=SLR_Function(model,seaLevel,sourceSink,x,heq0,heq,date,Vold,dxSLR,SLR0,SLROld)
%% Function or program name: SLR_Function
% SLR_FUNCTION G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       model — input argument (description/units to be verified)
%       seaLevel — input argument (description/units to be verified)
%       sourceSink — input argument (description/units to be verified)
%       x — input argument (description/units to be verified)
%       heq0 — input argument (description/units to be verified)
%       heq — input argument (description/units to be verified)
%       date — input argument (description/units to be verified)
%       Vold — input argument (description/units to be verified)
%       dxSLR — input argument (description/units to be verified)
%       SLR0 — input argument (description/units to be verified)
%       SLROld — input argument (description/units to be verified)
%
%   Outputs:
%       heq — output argument (description/units to be verified)
%       Vold — output argument (description/units to be verified)
%       model — output argument (description/units to be verified)
%       dxSLR — output argument (description/units to be verified)
%       SLR — output argument (description/units to be verified)
%       SLROld — output argument (description/units to be verified)
%       flagUpdate — output argument (description/units to be verified)
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

%SLR_FUNCTION G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.

% heq0 is the original equilibrium profile at t0
% heq is the evolved profile over some subsequent sea level rise.
% seaLevel is a structure containing sea level dates and values
% sourceSink is a structure containing volume /m coast of sediment
% additions or removal.
% Vold - previous sediment loss or addition m^3/m


% Get current values for SLR and FLUX
   SLR=interp1(seaLevel.dates, seaLevel.z,date,'linear','extrap')-SLR0;
   Vopt=interp1(sourceSink.dates, sourceSink.V,date,'linear','extrap');
   model.SLR=SLR;
   model.V=Vopt;
   flagUpdate=0;

% Only do this if there is a significant change
if abs(Vopt-Vold)>1 || abs(SLR-SLROld)>0.001
    disp('Computing SL adjustment')
     flagUpdate=1;
     heqOld=heq;
     SLROld=SLR; % Record old value to assess if change is large enough to warrent a sea level adjustment (NB total sea level rise)

% Compute translation
 [zout,dxSLR,q0,q,phi,dz]=computeTranslation(dxSLR,x,-heq,-heq0,SLR,model,Vopt-Vold);
 heq=-zout; % Key parameter - revised equilibrium parameter
 Vold=Vopt;% Record old flux

else
 %SLROld=SLR;
 %dxSLR=0;
end

return
end % SLR_Function

%% FUNCTIONS BELOW HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Shoreline translation subroutines (SLR and LS sediment fluxes)...
function [zout,dxSLR,q0,q,phi,dz]=computeTranslation(dxSLR,x,z,z0,SLR,model,q0opt)
%COMPUTETRANSLATION G_ForCE support routine; see function signature for inputs and outputs.

 options = optimset('TolX',1e-8');
 fun= @(dxSLR)transportCost(dxSLR,x,z,z0,SLR,-model.DOC,model.berm(end),q0opt,model);
 %[dxSLR,fval,exitflag,output] = fminsearch(fun,dxSLR,options);
 [dxSLR,fval,exitflag,output] = fminbnd(fun,-300,300,options);

 [q0,q,phi]=transportCost(dxSLR,x,z,z0,SLR,-model.DOC,max(model.berm),q0opt,model);

% Compute bed perturbation from volume change
 dz=gradient(q,x)';

% Avalanche if required
if model.avelancheOn==1
  jj=isnan(dz)==1;dz(jj)=0;
   zout=z+dz;
 [zout,iterations]=avelanche(x, zout,model);
else
 zout=z+dz;
end

% find NaNs
jj=find(isnan(zout)==1);
zout(jj)=z(jj);

end

 function [q0,q,phi]=transportCost(dxSLR,x,z,z0,SLR,DoC,Berm,q0opt,model)
%TRANSPORTCOST G_ForCE support routine; see function signature for inputs and outputs.

%% Make sure profile is a column vector
 x=x(:);z=z(:);z0=z0(:);

%% Define x-resolution and new x-scale
 xRes=model.transRes; % [m]
 xNew=min(x):xRes:max(x);

%% First interpolate profile at a finer resolution
 zNew=interp1(x,z0,xNew);
 zNew1=interp1(x,z,xNew);

%% Define z-offsets
 zPrime=zeros(size(xNew)); dz=zPrime;
 zPrime=zNew+SLR;

%% Compute spatially lagged profile at dx
 xPrime=xNew+dxSLR;
 zLag=interp1(xPrime,zPrime,xNew,'linear');

 %% Find active zone on profile
 jmin=find(zLag>=DoC+SLR,1,'first');
 jmax=find(zLag>=min(Berm+SLR,max(zLag)),1,'first');

 % Section added to develop a berm when sea level falls
 if SLR<0 || q0opt<0
  zLag(jmax:end)=Berm+SLR;
  jmax=find(zNew1>=Berm+SLR,1,'first');
  jmax=min(jmax,length(xNew));
  jmin=find(zNew1>=DoC+SLR,1,'first');
 end

  if model.BW==1
   kk=find(x<model.BWx);
   jmax=kk(end)-1;
  end

 jj=jmin:jmax;

 %% Compute Source function dz
 dz(jj)=zNew1(jj)-zLag(jj);
 dz=dz';

%% Call Poisson function to solve the Poisson Eq.
 % Define boundary type and values:
 % 0 - closed
 % 1 - free
 boundaryType=[1 0]; % 0 = closed (Neumann), 1 = free (Dirichlet)
 boundaryValues=[0 0]; % NB Specify a gradient for Neumann and a constant for Dirichlet
 [phi]=Poisson(xRes,zeros(size(dz)),dz,boundaryType,boundaryValues);

%% Compute sediment flux
 qTemp=gradient(phi,xRes);
 net=abs(qTemp(end)-qTemp(1)-q0opt); % Flux at boundary to be mimnimised to conserve volume
 gross=abs(mean(qTemp)+q0opt);
 if model.optMethod==0
     q0=net;
 else
     q0=gross;
 end

 %% Re-interpolate on original scale
  q=interp1(xNew,qTemp,x);
  phi=interp1(xNew,phi,x);

 end

%% Solve Poisson Eq.
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

%% Find solution via matrix inversion
function [y]=implicit(a,b,c,S)
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

%%
%% Avelanching algorithm
% Author Mark Davidson
function [zout,iterations]=avelanche(x,zin,model)
%AVELANCHE G_ForCE support routine; see function signature for inputs and outputs.

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
%  disp('************* Avelanche **************')

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
 boundaryType=[1 1]; % 0 = closed (Neumann), 1 = free (Dirichlet)
 boundaryValues=[zout(istart-1) zout(ifin+1)]; % NB Specify a gradient for Neumann and a constant for Dirichlet
%  boundaryType=[0 0]; % 0 = closed (Neumann), 1 = free (Dirichlet)
%  boundaryValues=[profGrad(istart) profGrad(ifin)]; % NB Specify a gradient for Neumann and a constant for Dirichlet
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
