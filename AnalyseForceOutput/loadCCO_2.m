function [chainage,z,dateNumber,RegID,OSGBx,OSGBy]=loadCCO_2
%% Function or program name: loadCCO_2
% Get input filename
%
%   Inputs:
%       None (script; see original code for workspace/file dependencies)
%
%   Outputs:
%       chainage — output argument (description/units to be verified)
%       z — output argument (description/units to be verified)
%       dateNumber — output argument (description/units to be verified)
%       RegID — output argument (description/units to be verified)
%       OSGBx — output argument (description/units to be verified)
%       OSGBy — output argument (description/units to be verified)
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

%% Get input filename
 disp('load CCO refeernce file')
 [FileName,PathName,FilterIndex] = ...
   uigetfile({'*.txt'},'Select beach profile data', 'MultiSelect','off');
 if FileName==0;RegID=0;chainage=0;dateNumber=0;z=0;return;end

%% Find datestring in filename
 istart=find(FileName=='_')+1;
 iend=istart+7;

%% Load file
 inputFile=[PathName FileName];% File and path
 [easting,northing,elevation,chainage,FC,profile,Reg_ID]...
   = textread(inputFile, ...
   '%f %f %f %f %s %s %s', -1,...
   'delimiter','\t','emptyvalue', NaN,'headerlines',1);

%% Look for unique profiles in the list
 uniqueProfiles=char(unique(Reg_ID));
 [n,m]=size(uniqueProfiles);
 disp(['Pick a profile number to analyse, 1 to ',num2str(n)])
% Display profile list
 for ii=1:n
  disp([num2str(ii) ') ',uniqueProfiles(ii,:)])
 end
 pickIndex=input('Pick a number and press return.... ');
 if isempty(pickIndex)==1; return; end
 profile2Analyse=uniqueProfiles(pickIndex,:);

%% Get date from filename
 dateStr=FileName(istart:iend);
 formatIn = 'yyyymmdd';
 dt=char(dateStr);
 dateNumber=datenum(dt,formatIn);
 disp(['Processing file recoded on the ',datestr(dateNumber)])

% Reg ID string
    S=char(Reg_ID);

%% Loop through all the Reg_IDs looking for the profile you want
    for ii=1:length(S);index(ii)=strcmp(S(ii,:),profile2Analyse);end
        jj=find(index==1);
% Check that the profile exists in the file
    if isempty(jj)==1; % Profile missing in file
        disp(['Profile missing in survey recorded on ',datestr(dateNumber)])
        qualityIndex=0;
    else % if exists then plot and extract contour
        OSGBx=easting(jj);
        OSGBy=northing(jj);
        xVar=chainage(jj);
        zVar=elevation(jj);
        [chainage,sortIndex]=sort(xVar);
        z=zVar(sortIndex);

        qualityIndex=1;
    end
        clear jj index

        RegID=profile2Analyse;
    end
