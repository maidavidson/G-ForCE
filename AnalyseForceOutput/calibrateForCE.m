%% Function or program name: calibrateForCE
% Programe to compare model output with field data
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
%% Programe to compare model output with field data
% Needs stats toolbox (although could be avoided)

%% initialise
%  clear
%  close all
 clc

%% Define calibration period
 calStart=[ ];%datenum(2008,1,1);
 calEnd =[ ];
 valStart=[ ]';
 valEnd=[ ];

%% Load field data and modeloutput
 disp('Please load field data file containing the structure survey')
 disp('survey.date')
 disp('survey.shoreline')
 disp('survey.volume')
 [fileName,filePath]=uigetfile('*.mat');
 if fileName==0;return;end
 filePathName=fullfile(filePath,fileName);
 load(filePathName)

 %% Load model output file
 disp(' ')
 disp('Now load the ForCE model output file for comparison')
 [fileName,filePath]=uigetfile('*.mat');
 filePathName2=fullfile(filePath,fileName);
 load(filePathName2)

%% Iterpolate model on same time-scale as field data
 modelShoreline=interp1(model.dateArray,model.contour,survey.date);
 modelVolume=interp1(model.dateArray,model.V,survey.date);

%% Find calibration values
 if isempty(calStart); calStart=min(survey.date);end
 if isempty(calEnd); calEnd=0.5*(max(survey.date)+min(survey.date));end
 jj=find(survey.date>=calStart & survey.date<=calEnd);

 %% Find validation values
 if isempty(valStart); valStart=calEnd;end
 if isempty(valEnd); valEnd=max(survey.date);end
 kk=find(survey.date>=valStart & survey.date<=valEnd);

%% Calibration regression between data and model
 [Rsc,ms,bs] = regression(modelShoreline(jj),survey.shoreline(jj));
 [Rvc,mv,bv] = regression(modelVolume(jj),survey.volume(jj));
 [Rsv] = regression(modelShoreline(kk),survey.shoreline(kk));
 [Rvv] = regression(modelVolume(kk),survey.volume(kk));
 [Rs,mts,bts] = regression(modelShoreline,survey.shoreline);
 [Rv,mtv,btv] = regression(modelVolume,survey.volume);

%% Legend text
legend1=['Model'];
legend2=['Calibration Data: r = ',num2str(Rsc)];
legend3=['Validation Data: r = ',num2str(Rsv)];

legendV1=['Model'];
legendV2=['Calibration Data: r = ',num2str(Rvc)];
legendV3=['Validation Data: r = ',num2str(Rvv)];

%% Calibrate model
 shorelineCalModel=ms*model.contour+bs;
 volumeCalModel=mv*model.V+bv;

 %%
 figure(2);
 clf
 set(gcf,'units','normalized','outerposition',[0 0 1 1]);
 subplot(2,1,1)
 plot(model.dateArray,shorelineCalModel,'Color',[0. 0. 0.],'LineWidth',2,'LineStyle','-','Marker','none')
 hold on
 plot(survey.date(jj),survey.shoreline(jj),'o','MarkerFaceColor','k','MarkerEdgeColor','k')
 plot(survey.date(kk),survey.shoreline(kk),'s','MarkerFaceColor','w','MarkerEdgeColor','k')
 grid on
 legend(legend1, legend2, legend3)
 ylabel('Displacement [m]')
 title(['Model calibration: multiply model.k1 by ', num2str(mts),', r = ',num2str(Rs)])
 datetick
 xlim([min(survey.date) max(survey.date)])
 set(gca,'FontSize',20,'LineWidth',3);
 %
 subplot(2,1,2)
 plot(model.dateArray,volumeCalModel,'Color',[0. 0. 0.],'LineWidth',2,'LineStyle','-','Marker','none')
 hold on
 plot(survey.date(jj),survey.volume(jj),'o','MarkerFaceColor','k','MarkerEdgeColor','k')
 plot(survey.date(kk),survey.volume(kk),'s','MarkerFaceColor','w','MarkerEdgeColor','k')
 grid on
 legend(legendV1, legendV2, legendV3)
 title(['Model calibration: multiply model.k1 by ', num2str(mtv),', r = ',num2str(Rv)])
 xlabel('Date')
 ylabel('Volume [m^3/m coast]')
 datetick
 xlim([min(survey.date) max(survey.date)])
 set(gca,'FontSize',20,'LineWidth',3);


 return

 %%  NMSE
  Ev=rms(calModel(kk)-survey.x(kk)).^2./std(survey.x).^2;
  Ec=rms(calModel(jj)-survey.x(jj)).^2./std(survey.x).^2;

 %% Combined score
CSS=(Rv(1,2).^2+Rc.^2)/(Ev+Ec);

 disp(['Model Cal. r-value = ',num2str(Rc)])
 disp(['Model Val. r-value = ',num2str(Rv(1,2))])
 disp(['Model Cal. nrms-value = ',num2str(Ec)])
 disp(['Model Val. nrms-value = ',num2str(Ev)])
 disp(['Model CSS-value = ',num2str(CSS)])
 disp(['Model m-value = ',num2str(m)])
 disp(['Model b-value = ',num2str(b)])
