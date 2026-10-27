function [tArray,zArray, model, profile]=getModelOutput(fname,dateStart,dateFinish)
%% Function or program name: getModelOutput
% Get ForCE output datafile
%
%   Inputs:
%       fname — input argument (description/units to be verified)
%       dateStart — input argument (description/units to be verified)
%       dateFinish — input argument (description/units to be verified)
%
%   Outputs:
%       tArray — output argument (description/units to be verified)
%       zArray — output argument (description/units to be verified)
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

%% Get ForCE output datafile
if exist(fname,'file')==2
    load(fname)
 else
     disp('Output file not found');return;
end

%% Find file output file
 S=whos('-file', fname);
 SS=char(S.name);

%% load in profiles and dates for all years stored in the file
 kk=1;
 for ii =1:length(S)
   index(ii)=startsWith(SS(ii,:),'Profiles');
   if index(ii)==1
       SS2=strtrim(SS(ii,:));
       variableName(kk,:)=SS(ii,:);
      % mm=length(SS2(1,:));
       dateVal=SS2(end-3:end);
       numDate(kk)=str2num(dateVal);
       kk=kk+1;
   end
end

%% Sort in date order
 numDate=sort(numDate);

%% Temp dnum
 dnumTemp=datenum(numDate,1,1);
 jj=find(dnumTemp>dateStart & dnumTemp<dateFinish);

%% Find number of years of data
 nYears=length(jj);
 if nYears==0
     warning('getModelOutput:NoYears', 'No saved annual outputs fall within the requested date range.');
     zArray=[];
     tArray=[];
     load(fname,'profile','model');
     return
 end

%% Initialise elevation and time arrays
 zArray=[];
 tArray=[];

%% Load the model structure from the file
 load(fname,'profile','model');

 %% Define time-resolution of plot
  nPts=1; % Plot frequency (plot step = nPts)
  Res=nPts*model.dateIncrement;

%% Loop through all years extracting data into arrays
 for ii=jj(:)'
  variable1=['Profiles_',num2str(numDate(ii))];
  variable2=['Dates_',num2str(numDate(ii))];
  load(fname,variable1,variable2);
  eval(['zArray=[zArray,',variable1,'(:,1:Res:end)];']);
  eval(['tArray=[tArray,',variable2,'(1:Res:end)];']);
  eval(['clear ',variable1,' ',variable2])
 end

 end
