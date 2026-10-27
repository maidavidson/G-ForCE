%% Function or program name: prepare_CCO_4_ForCE
% Programe to prepare CCO profile for ForCE
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
%% Programe to prepare CCO profile for ForCE
% This makes the appropriate structure and extrapolates the profile
% to beyound the depth of closure (zExtrap) using the average
% gradient of the lower beach profile below zMax.

%% Initialise
 %clear all
 %close all
 clc

%% Load CCO topo datafile
 [chainage,z,dateNumber,RegID]=loadCCO;
 survey.x=chainage; survey.x= survey.x(:);
 survey.z=z; survey.z= survey.z(:);
 survey.date= dateNumber;
 survey.RegID= RegID;

%% Make survey filename
 opFile=['surveyData_',datestr(dateNumber)];
 uisave({'survey'}, opFile)

 % Plot data
    %% Plot Profiles
    fig1=figure;
    axes1 = axes('Parent',fig1);
    hold(axes1,'all');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);

    plot(survey.x,survey.z)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    legend(datestr(dateNumber))
    grid on
    title(['Profile Reg ID = '  survey.RegID])
    xlabel('Chainage (m)')
    ylabel('Beach profile height (m ODN)')
    h=gca;
    set(h,'FontSize',16);
