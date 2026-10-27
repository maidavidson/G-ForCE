function [zPrime,q,D,DSand,model,maxRunup_x,maxWaterDepth,Gradient,swashIndex,Ir,overwash]=...
    ForCEfunction3(x,heq,zPrime,model,constant,wave,tide,seaLevel,dates,jRock,dxSLR,SLR,nrep)
%% Function or program name: ForCEfunction3
% Compute instantaneous dissipation
%
%   Inputs:
%       x — input argument (description/units to be verified)
%       heq — input argument (description/units to be verified)
%       zPrime — input argument (description/units to be verified)
%       model — input argument (description/units to be verified)
%       constant — input argument (description/units to be verified)
%       wave — input argument (description/units to be verified)
%       tide — input argument (description/units to be verified)
%       seaLevel — input argument (description/units to be verified)
%       dates — input argument (description/units to be verified)
%       jRock — input argument (description/units to be verified)
%       dxSLR — input argument (description/units to be verified)
%       SLR — input argument (description/units to be verified)
%       nrep — input argument (description/units to be verified)
%
%   Outputs:
%       zPrime — output argument (description/units to be verified)
%       q — output argument (description/units to be verified)
%       D — output argument (description/units to be verified)
%       DSand — output argument (description/units to be verified)
%       model — output argument (description/units to be verified)
%       maxRunup_x — output argument (description/units to be verified)
%       maxWaterDepth — output argument (description/units to be verified)
%       Gradient — output argument (description/units to be verified)
%       swashIndex — output argument (description/units to be verified)
%       Ir — output argument (description/units to be verified)
%       overwash — output argument (description/units to be verified)
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

% Compute instantaneous dissipation
if isfield(model,'MorphFB')
    hpass=heq-model.MorphFB*zPrime;
else
    hpass=heq;
end

[D,DSand,~,overwash,maxRunup_x,maxWaterDepth,Gradient,swashIndex,Ir]=...
    findMeanDissipation_3(x,hpass,model,wave,tide,seaLevel,dates,model.dt);
Dmax=max(abs(D),abs(DSand));

% Find layer positions based on pre-defined depths
if model.layers == 2
    % The gravel boundary must be higher than the mixed-zone boundary.
    % Reject invalid ordering or boundaries outside the computational profile.
    if model.gravel <= model.mixed
        error('G_ForCE:InvalidLayerBoundaries', ...
            'model.gravel (%.3g) must exceed model.mixed (%.3g).', ...
            model.gravel, model.mixed);
    end
    surfaceElevation = -heq + zPrime;
    igravel = find(surfaceElevation > model.gravel + SLR, 1, 'first');
    imixed = find(surfaceElevation > model.mixed + SLR, 1, 'first');
    if isempty(igravel) || isempty(imixed)
        error('G_ForCE:LayerBoundaryOutsideDomain', ...
            ['A two-layer boundary lies outside the computational domain. ' ...
             'Check model.gravel, model.mixed, SLR and the profile extent.']);
    end
    if imixed > igravel
        error('G_ForCE:LayerBoundaryOrder', ...
            'Mixed-zone index must not exceed gravel-zone index.');
    end
end

% Find depth of closure
    idoc=find(-(heq-zPrime)>-(model.DOC+SLR),1,'first');
    if isempty(idoc); idoc=1;end

% Translate mean dissipation map
    Deq=interp1(x(:)+dxSLR,model.Deq,x(:),'linear','extrap');
    Deq2=interp1(x(:)+dxSLR,model.Deq2,x(:),'linear','extrap');
    q1=zeros(size(Deq));
    q2=q1;
    zPrime1=zeros(size(Deq));
    zPrime2=zeros(size(Deq));

    %%%%%%%%% generate window functions
    if model.layers==1

        w1=ones(size(q1)); % Initialise window for gravel transport
        w2=w1.*0;
        w1(1:idoc)=0;

    else

        w1=ones(size(q1)); % Initialise window for gravel transport
        w2=w1; % Initialise window for sand transport
        indexMixed=imixed:igravel; % find mixed beach index
        nTaper=length(indexMixed);
        [w2(indexMixed),w1(indexMixed)]=cosexpdecay(nTaper);
        w2(1:idoc-1)=0; % Nullify area offshore of DOC
        w2(igravel+1:end)=0;
        w1(1:imixed-1)=0;

    end
    %%%%%%%%%% end window definition


% Compute cross-shore sediment flux in all layers
    if Ir>=model.IrThreshold || model.berm(end)==model.tideThreshold
        q1 =  w1.* model.k1(1).*constant.k1.*(D) ; % gravel
    else
        q1 =  w1.*model.k1(1).*constant.k1.*(Deq-D) ; % gravel
    end

    if ~isempty(model.k1(2))
        q2 =  w2.*model.k1(2).*constant.k1.*(Deq2-DSand); % sand
    end


% adjust transport if overwash parameter is not still zero...
if overwash ~=0
    q1 = abs(w1 .*model.k1(1).*constant.k1.*D) ;
    model.nOverwash=model.nOverwash+1;
end


% Set q to near zero where there's rock/non-erodible bed
if ~isempty(jRock)
    q1(jRock)=0.0001*q1(jRock);
    q2(jRock)=0.0001*q2(jRock);
end


% Zero inactive q calculations where Deq is finite and D=0.
    % Gravel
    jj=find(abs(D)==0);
    q1(jj)=0;
    if model.layers==2
        q1(1:imixed)=0;
    end
    % Sand
    jj=find(abs(DSand)==0);
    q2(jj)=0;
    q2(1:idoc-1)=0;
    if model.layers==2
        q2(igravel+1:end)=0;
    end
% total flux
q= q1+q2;

% Compute bed perturbation (Exner equ.)
if max(abs(q))~=0
   zPrime1 = -gradient(q1,x(:))'*model.dt*constant.days2seconds; % Gravel
   zPrime2 = -gradient(q2,x(:))'*model.dt*constant.days2seconds; % Sand
end


% Combine perturbation components
zPrime=zPrime(:) + zPrime1(:) + zPrime2(:);zPrime=zPrime(:)';

jj=find(Dmax>0); % Find non-zero dissipation

if max(abs(zPrime))~=0 & ~isempty(jj)% do diffusion

% Define diffusion function
 Diss=zeros(size(Dmax)); % Initialize diffusion coeff.
 Dm=mean(abs(Dmax(jj))); % Compute max value
 Dm=max(Dm,1); % make sure that Dm is non-zero and +ve
 endi=length(x);
 if isempty(swashIndex) | isempty(jj)
     indexNeeded=idoc:idoc+1;
 else
     indexNeeded=idoc:min(max(swashIndex(end)+1,jj(end)),endi);
 end
 Diss(indexNeeded)=model.k2(1) * constant.k1 * Dm ./ Gradient; % compute dissipation coefficient for active zone
 Diss=Diss';

% Diffusion algorithm
if length(jj)>5

  % Diffusion boundary conditions
  boundaryType  =[0 0]; % 0=Zero Flux, 1 = pinned
  boundaryValues=[0 0]; % NB Specify a gradient for Neumann and a constant for Dirichlet

[zPrime]=crankNicholson...
        (zPrime,Diss,model.dt*constant.days2seconds,...
         model.dx,model.fact,boundaryType,boundaryValues);
end

end % End do diffusion logic


end % End subroutine
