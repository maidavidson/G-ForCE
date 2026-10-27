%% Function or program name: plotModelProfileEvolution
% Get output file name
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
% Get output file name
 %% Get setup file
 disp('Please load the model output file')
 [file,path] = uigetfile('*.mat');

%% Now load output data and plot
fname=fullfile(path,file);
dateStart=datenum(1980,1,1);
dateFinish=datenum(2100,1,1);
[tArray,zArray, model, profile]=getModelOutput(fname,dateStart,dateFinish);

figure
subplot(2,1,1)
plot(model.x,zArray,'k','LineWidth',3);
hold on
plot(model.x,-profile.heq,'w.','LineWidth',3);
title(['Profile ID: ',profile.ID], 'Interpreter', 'none')
ylim([-model.DOC,1.1*max(-profile.h0)])
jx=find(profile.heq<=model.DOC,1,'first');
xlim([model.x(jx) model.x(end)])
grid on
box on
set(gca,'LineWidth',3,'FontSize',24)

subplot(2,1,2)
zplot=zArray'+profile.h0; % plot change rel. to starting profile
jj=find(abs(zplot)<0.02);
zplot(jj)=0;
cHandle=contourf(model.x,tArray,zplot);
hold on
colormap(hot)
xlim([model.x(jx) model.x(end)])
datetick('y')
c=colorbar;
c.Location='northoutside';
set(gca,'LineWidth',3,'FontSize',24)
