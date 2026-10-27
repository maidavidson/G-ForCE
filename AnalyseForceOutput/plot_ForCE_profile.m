%% Function or program name: plot_ForCE_profile
% Plot model profile
%
%   Inputs:
%       None (script; see original code for workspace/file dependencies)
%
%   Outputs:
%       None (no formal return arguments)
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
%% Plot model profile
% program selects a profile from the input data file and saves it a
% structure called profile

%% Initialise
 %clear all
 %close all
 clc

%% Input profile date needed
 disp('Please pick the profile date needed. Make sure that it is in the modelled range!')
 h = uicontrol('style', 'edit', 'position', [10 10 100 20]);
 uicalendar('DestinationUI', {h, 'String'});

%% Now wait for the string to be updated
 waitfor(h,'String');
 dateRequired = datenum(get(h,'String'));

%% Load profile from chosen ForCE o/p datafile
 disp('Load ForCE output datafile')
 [filename,filepath] = uigetfile('*.mat','Tite','Load ForCE output datafile');
 if filename==0;return;end
 [x,z0,ze,zPrime,model,profile]=getModelProfile(dateRequired,filepath,filename);
 [dum1,dum2,dum3,zPrime0,dum4,dum5]=getModelProfile(profile.t0,filepath,filename);

%% Plot results
 modelProfile.x=x;
 modelProfile.z=ze+zPrime;
 modelProfile.z0=ze+zPrime0;
 modelProfile.date=dateRequired;
 modelProfile.dT=dateRequired-profile.t0;

 plot(modelProfile.x,modelProfile.z,'k','LineWidth',2)
 hold on
 plot(modelProfile.x,modelProfile.z0,'b:','LineWidth',2)
 grid on
 box on
 xlabel('Cross-shore distance (m)')
 ylabel('Beach profile height [m]')
 set(gca,'FontSize',20,'LineWidth',3);
 title(['Beach profile for ',model.site,'. {\Delta}T = ',num2str(modelProfile.dT),' days.'])
 legend(['Prediction: ',datestr(modelProfile.date)],['Initial Profile: ',datestr(profile.t0)])

%% Save selected model output
 disp('Please enter a name for the output profile file')
 opProfile=['opProfile_',datestr(dateRequired)];
 uisave({'modelProfile','model'},opProfile)
