%% Function to extract a model profile from the output datafile
% Input
% dateRequired [datenum format]
% filepath - path to output datafile
% filename - output datafilename
% Output
% x - cross-shore distance [m] increasing seawards
% z0 - initial profile (used to initialise model at t0), zero at MSL
% ze - equilibrium profile
% xPrime - profile pertibation rel. to ze
% NB h(x)=-(ze+zPrime)
% model - structure containing all model parameters

function [x,z0,ze,zPrime,model,profile]=getModelProfile(dateRequired,filepath,filename)
%% Function or program name: getModelProfile
% Function to extract a model profile from the output datafile
%
%   Inputs:
%       dateRequired — input argument (description/units to be verified)
%       filepath — input argument (description/units to be verified)
%       filename — input argument (description/units to be verified)
%
%   Outputs:
%       x — output argument (description/units to be verified)
%       z0 — output argument (description/units to be verified)
%       ze — output argument (description/units to be verified)
%       zPrime — output argument (description/units to be verified)
%       model — output argument (description/units to be verified)
%       profile — output argument (description/units to be verified)
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


 [yr,mm,dd,hr,minute,sec]=datevec(dateRequired);

%% Initialise elevation and time arrays
 zArray=[];
 tArray=[];

%% Find file output file
 S=fullfile(filepath,filename);

%% Load in the model dates for the relevant year
 variable2=['Dates_',num2str(yr)];
 load( S,variable2);
 load( S,'model');
 eval(['tArray=[tArray,',variable2,'(1:end)];']);
 minimise=min(abs(tArray-dateRequired));
 jj=find( abs(tArray-dateRequired) == minimise);

%% Load in the relevant years profiles
 variable1=['Profiles_',num2str(yr)];
 load( S,variable1);
 eval(['zArray=[zArray,',variable1,'(:,jj)];']);

 eval(['clear ',variable1,' ',variable2])

 %% Load the model structure from the file
 load( S,'profile','model');
 x=-profile.x;
 zPrime=zArray';
 z0=-profile.h0';
 ze=-profile.heq;
end
