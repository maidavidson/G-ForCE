function [berm]=bermCalculator(wave,tide,model,SLR)
%% Function or program name: bermCalculator
% BERMCALCULATOR G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       wave — input argument (description/units to be verified)
%       tide — input argument (description/units to be verified)
%       model — input argument (description/units to be verified)
%       SLR — input argument (description/units to be verified)
%
%   Outputs:
%       berm — output argument (description/units to be verified)
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

%BERMCALCULATOR G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
disp('Computing Berm Height')
 minT=max(min(wave.dates),min(tide.dates));
 maxT=min(max(wave.dates),max(tide.dates));
 dt=1/(24*4);
 timeScale=minT:dt:maxT;
 tideE  =interp1(tide.dates,tide.z           ,timeScale);
 tideE=tideE+SLR;
 H=interp1(wave.dates,wave.H,timeScale);
 Tp=interp1(wave.dates,wave.Tp,timeScale);
 [R2]=runupSimple(H,0.17);
 R=tideE+R2;
 berm=max(R);
 disp(['Finished Computing Berm Height : ',num2str(berm), 'm'])
end
