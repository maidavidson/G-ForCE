function [xout,zout]=extrapolateSurvey(chainage,z,zExtrap,zMax,modelRes)
%% Function or program name: extrapolateSurvey
% EXTRAPOLATESURVEY G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       chainage — input argument (description/units to be verified)
%       z — input argument (description/units to be verified)
%       zExtrap — input argument (description/units to be verified)
%       zMax — input argument (description/units to be verified)
%       modelRes — input argument (description/units to be verified)
%
%   Outputs:
%       xout — output argument (description/units to be verified)
%       zout — output argument (description/units to be verified)
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

%EXTRAPOLATESURVEY G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% chanage - distance down the beach in offshore direction
% z - beach level usually rel. MSL
% zExtrap - max profile z needed should be < DOC e.g. = -30m
% zMax - z-level below which the gradient is computed
% modelRes - Model resolution e.g. =5;

% Number of survey points
 N=length(z);

% Compute gradient
 grad=gradient(z,chainage);

% Find active beach gradient
 jj=find(z>=zMax,1,'last');

% Mean gradient
 meanGrad=mean(grad(jj:end));

% distance to max depth
 depthChange=zExtrap-min(z);
 chainageChange=depthChange./meanGrad;

% Add deep point to profile and interpolate onto regular grid
 chainage(N+1)=chainage(N)+ chainageChange;
 z(N+1)= zExtrap;
 xout=min(chainage):modelRes:max(chainage); xout= xout(:);
 zout=interp1(chainage,z,xout,'linear'); zout= zout(:);
end
