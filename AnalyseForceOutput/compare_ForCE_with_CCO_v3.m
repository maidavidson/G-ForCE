%% Function or program name: compare_ForCE_with_CCO_v3
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
% Load a reference survey, then a subsequent survey
% Finally load a model output file to compare with the subsequent survey
% file

%% Clear workspace
 clear all
 close all
 clc

%% Instructions
disp('You will be asked to load 2 files')
disp('1) Spin-up survey file')
disp('2) The ForCE output data file')
disp('')

% Load ref profile
% disp('load the first CCO reference profile (t0)')
% [filename,filepath] = uigetfile;
% load (strcat(filepath,filename))

load('/Users/mdavidson/Documents/ForCE_G/inputData/Borth/Borthprofile13LIDARspinUp.mat')

Nprofile=length(survey.date);
chainage1=survey.x;
z1=survey.zArray(:,1)';

% % Load comparison profile
% disp('Load the later CCO comparison profile')
% [chainage2,z2,dateNumber2,RegID2]=loadCCO;


LineType{1}='--';
LineType{2}=':';
LineType{3}='-.';


%% Get model data file
disp('Load FoRCE output file for comparison')
 [filename,filepath] = uigetfile('MultiSelect','on');

for j=1:length(filename)

% get first profile and interpolate on the same chainage as survey
 if iscell(filename)
     [x,z0,ze,zPrime1,model,profile]=getModelProfile(survey.date(1),filepath,cell2mat(filename(j)));
 else
     [x,z0,ze,zPrime1,model,profile]=getModelProfile(survey.date(1),filepath,filename);

 end

 for i = 2:Nprofile

 if iscell(filename)
     [x,z0,ze,zPrime2,model,profile]=getModelProfile(survey.date(i),filepath,cell2mat(filename(j)));
 else
     [x,z0,ze,zPrime2,model,profile]=getModelProfile(survey.date(i),filepath,filename);

 end

legendStr{1,i-1}=['Measured initial profile: ',datestr(survey.date(1))];
legendStr{2,i-1}=['Measured subsequent profile: ',datestr(survey.date(i))];

 z2=survey.zArray(:,i)';
%% Compute relative change
 delta_z=interp1(x,(zPrime2-zPrime1),chainage1);
 zm=survey.zArray(:,i)';
 predicted_z=z1+delta_z;
 % predicted_z=interp1(x,zPrime2,chainage1);

istart=1;
maxPlot=4;
maxPlot0=4;
%% Plot Profiles
    if i-1>maxPlot
        istart=istart+1;
        maxPlot=maxPlot+maxPlot0;
    end
    figure(istart)
    subplot(2,2,i-1);
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    if j==1
        plot(chainage1,z1,'b','LineWidth',3)
        hold on
        plot(chainage1,z2,'r','LineWidth',3)

    end
    % jj=find(predicted_z>=model.berm,1,'last');
    % predicted_z(1:jj)=z1(1:jj);

    plot(chainage1,predicted_z,'k','LineWidth',2,'LineStyle',LineType{j})
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    legend(['Measured initial profile: ',datestr(survey.date(1))],['Measured subsequent profile: ',datestr(survey.date(i))],'','','');

    grid on
     box on
    xlabel('Cross-shore distance (m)')
    ylabel('Beach profile height [m]')
    h=gca;
    set(h,'FontSize',16,'LineWidth',3);

% Skills
Berm=5;
jj=find(~isnan(zm)  & z1<Berm);
ME=0.1;
[BSS(i-1)]=brierss(zm(jj),predicted_z(jj),z1(jj),ME);
[rmsError(i-1)]=rms(predicted_z(jj)-zm(jj));


legendStr{2+j,i-1}=['k_1 = ',num2str(model.k1(1)),' BSS = ', num2str(BSS(i-1)),' R.M.S. error = ', num2str(rmsError(i-1)),'m'];
title(['Model prediction: ',datestr(survey.date(i)),', \DeltaT = ',num2str(survey.date(i)-survey.date(1)),' days'])

end

end % End files

subplot(2,2,1)
legend(legendStr{:,1})
subplot(2,2,2)
legend(legendStr{:,2})
subplot(2,2,3)
legend(legendStr{:,3})
subplot(2,2,4)
legend(legendStr{:,4})


disp('R.M.S. error = ')
mean(rmsError)
disp('Brier Skills Score = ')
mean(BSS)

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
