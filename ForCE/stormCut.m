function  [heq]=stormCut(x,heq,tide,model)
%% Function or program name: stormCut
% STORMCUT G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       x — input argument (description/units to be verified)
%       heq — input argument (description/units to be verified)
%       tide — input argument (description/units to be verified)
%       model — input argument (description/units to be verified)
%
%   Outputs:
%       heq — output argument (description/units to be verified)
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

%STORMCUT G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.

% remember old heq
heqOld=heq;

%% Compute limits
MHW=max(tide.z);
MWL=mean(tide.z);
i1=find(-heq>=MWL,1,'first');
i2=find(-heq>=MHW,1,'first');
profileGradient=gradient(heq,x);
shorelineGradient= max(abs(profileGradient(i1:i2)));

%% Generate modified heq for region above HW to swash limit
 k1=find(-heq>=MHW,1,'first');
 k2=length(heqOld);
 gradheq=gradient(-heq,x);
 gradheq(k1:k2)=shorelineGradient;
 hNew=-cumtrapz(x,gradheq)+heq(1);
 heq(k1:end)=hNew(k1:end);
 heq=max(heq,heqOld);
 [zout,iterations]=avelanche(x,-heq,model);
 heq=-zout;
end
