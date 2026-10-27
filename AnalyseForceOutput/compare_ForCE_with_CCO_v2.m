%% Function or program name: compare_ForCE_with_CCO_v2
% 9/10/20
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

%% Written by Mark Davidson
%  9/10/20
%  Program reads CCO data files and ForCE model output for comparison
%  Need to download CCO topo files
% Load a reference survey, then a subsequent survey
% Finally load a model output file to compare with the subsequent survey
% file

%% Clear workspace
 clear all
 close all
 clc

%% Instructions
disp('You will be asked to load 3 files')
disp('1) CCO refernce file (NB The date should match the start of the model run)')
disp('2) CCO comparison file (recorded after the reference file)')
disp('3) The ForCE output data file')
disp('')

% Load ref profile
disp('load the first CCO reference profile (t0)')
[chainage1,z1,dateNumber1,RegID1]=loadCCO;
if RegID1==0;return;end

% Load comparison profile
disp('Load the later CCO comparison profile')
[chainage2,z2,dateNumber2,RegID2]=loadCCO;

%% Get model data file
disp('Load FoRCE output file for comparison')
 [filename,filepath] = uigetfile;
 [x,z0,ze,zPrime1,model,profile]=getModelProfile(dateNumber1,filepath,filename);
 [x,z0,ze,zPrime2,model,profile]=getModelProfile(dateNumber2,filepath,filename);

%% Compute relative change
 delta_z=interp1(x,(zPrime2-zPrime1),chainage1);
 zm=interp1(chainage2,z2,chainage1);
 predicted_z=z1+delta_z;
 % predicted_z=interp1(x,zPrime2,chainage1);

%% Plot Profiles
    figure(2);
    clf
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    plot(chainage1,z1,'b','LineWidth',3)
    hold on
    plot(chainage2,z2,'r','LineWidth',3)
    % jj=find(predicted_z>=model.berm,1,'last');
    % predicted_z(1:jj)=z1(1:jj);
    plot(chainage1,predicted_z,'k--','LineWidth',2)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    legend(['Measured initial profile: ',datestr(dateNumber1)],['Measured subsequent profile: ',datestr(dateNumber2)],...
        ['Model prediction: ',datestr(dateNumber2),', \DeltaT = ',num2str(dateNumber2-dateNumber1),' days'])%,...
    grid on
    title(['Profile Reg ID = ' RegID1, ', Prediction duration: ', num2str(dateNumber2-dateNumber1),' days'])
     box on
    xlabel('Cross-shore distance (m)')
    ylabel('Beach profile height [m]')
    h=gca;
    set(h,'FontSize',24,'LineWidth',3);

% Skills
Berm=4;
jj=find(~isnan(zm)  & z1<Berm);
disp('Brier Skills Score = ')
ME=0.1;
[BSS]=brierss(predicted_z(jj),zm(jj),z1(jj),ME)
disp('R.M.S. error = ')
[rmsError]=rms(predicted_z(jj)-zm(jj))


return
  %% Skills
function [BSS]=brierss(xm,x,xb,dx)

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% Brier Skils Function

% xm - measured data
% x - model data
% xb - benchmrk data
% dx - measurement error

BSS = 1 - ...
    ( ...
    mean( ( ( abs(x-xm) - dx ) ).^2) /...
    mean( (xb-xm).^2) ...
    );
end
