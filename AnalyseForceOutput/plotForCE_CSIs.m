%% Function or program name: plotForCE_CSIs
% Program to plot ForCE CSIs
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
%% Program to plot ForCE CSIs
% These include the tanslation of the shoreline contour
% and intertidal volume time-series
% The vertical elevation of the shorelne contour and volume
% calculations are defined in the setup file

%% Initialise
 %clear all
 %close all
 clc

%% Get ForCE output datafile
disp('Load ForCE output datafile')
 [filename,path] = uigetfile('*.mat','Title','Load ForCE output file');
 if filename==0;return;end

%% Define date range to be plotted
 dateStart=datenum(1980,1,1);
 dateFinish=datenum(2100,1,1);

%% Find file output file
 filePathName=fullfile(path,filename);

%% Load the model structure from the file
 load(fullfile(path,filename),'model','wave');


%% Plot
   figure(2);
   clf
   set(gcf,'units','normalized','outerposition',[0 0 1 1]);
  subplot(2,1,1)
  plot(wave.dates, wave.H,'k','LineWidth',2)
  datetick
  xlim([min(model.dateArray), max(model.dateArray)])
  grid on
  box
  ylabel('Wave height [m]')
  title(['ForCE CSIs for fieldsite: ',model.site])
  haxis=gca;
  set(haxis,'FontSize',16);

 subplot(2,1,2)
 % LH axis
  yyaxis left;
  plot(model.dateArray,model.contour-model.contour(1),'k','LineWidth',2')
  hold on
  plot(model.dateArray,model.maxShoreline-model.maxShoreline(1),'r:','LineWidth',2)
  xlim([min(model.dateArray), max(model.dateArray)])
  hold off
  ax=gca;
  ax.YColor='black';
  ylabel('Displacement [m]')
 % RH axis
  yyaxis right;
  plot(model.dateArray,model.V-model.V(1),'b','LineWidth',2)
  ylabel('Beach Perturbation [m]')
  ax=gca;
  ax.YColor='black';
  datetick('x','yyyy')
  xlim([min(model.dateArray), max(model.dateArray)])
  grid on
  xlabel('Date')
  ylabel('Volume [m^3/m coast]')
%   title(['Model contour: z = ', num2str(model.zContour),'m. ',...
%       'Model volume: z = ', num2str(model.zmin),...
%        ' to ', num2str(model.zmax), 'm. '...
%       'Last Update: ',  datestr(model.date),'.' ])
  legend('Shoreline Displacement','Maximum Shoreline Displacement','Volume','Location','northwest')
  haxis=gca;
  set(haxis,'FontSize',16);
