%% Function or program name: processCCO
% 1/3/2019
%
%   Inputs:
%       None (script; see original code for workspace/file dependencies)
%
%   Outputs:
%       Saved model results (see save calls in script)
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
%  1/3/2019
%  Program reads CCO data files, saves volume, shoreline and profiles
%  Need to download CCO topo files including baseline files
% Fix missing profile ID bug
% Fix volume calculation bug

%% Clear workspace
    clear
    close all
    clc

%% Initialise profile reg id array
    Reg_IDArray=[];

%% Get input filename
   disp('Please multi-select (hold shift key) all CCO topo files')
    [FileName,PathName,FilterIndex] = ...
        uigetfile({'*.txt'},'Select beach profile data', 'MultiSelect','on');

if iscell(FileName)==0; return;end

%% Determnine if one file or multiple files
    if ischar(FileName)
        NN=1;
    else
        NN=length(FileName);
    end

% Find datestring in filename
    istart=find(char(FileName(1))=='_')+1;
    iend=istart+7;

%% Loop through all the files selected
    for kk=1:NN % Loop through all files

% Load first file and determine what profiles are in them
    if ischar(FileName)
        fName=FileName;
        inputFile=[PathName FileName];
    else
        fName=char(FileName(1));
        inputFile=[PathName char(FileName(kk))] ;
    end

% Pick profile to analyse
    [eastin,northing,elevation,chainage,FC,profile,Reg_ID]...
        = textread(inputFile, ...
        '%f %f %f %f %s %s %s', -1,...
        'delimiter','\t','emptyvalue', NaN,'headerlines',1);
        Reg_IDArray=[Reg_IDArray; Reg_ID];
    end % end file loop

%% Look for unique profiles in the list
    uniqueProfiles=char(unique(Reg_IDArray));
    [n,m]=size(uniqueProfiles);
    disp(['Pick a profile number to analyse, 1 to ',num2str(n)])

 % Display profile list
    for ii=1:n
        disp([num2str(ii) ') ',uniqueProfiles(ii,:)])
    end
    pickIndex=input('Pick a number and press return.... ');
    profile2Analyse=uniqueProfiles(pickIndex,:);
    profile2Analyse=strtrim(profile2Analyse); % added 26/1/22 removes white space


%% Now loop through all files again extracting the profile requested
    nProfile=1; % No successful profile extractions.
    for mm=1:NN
        if ischar(FileName)
            fName=FileName;
            inputFile=[PathName FileName];
        else
            fName=char(FileName(mm));
            inputFile=[PathName char(FileName(mm))];
        end

% Get date from filename
        dateStr(mm,:)=fName(istart:iend);
        formatIn = 'yyyymmdd';
        dt=char(dateStr(mm,:));
        dateNumber(mm)=datenum(dt,formatIn);
        disp(['Processing file recoded on the ',datestr(dateNumber(mm))])

% Display file being loaded
        [easting,northing,elevation,chainage,FC,profile,Reg_ID]...
            = textread(inputFile, ...
            '%f %f %f %f %s %s %s', -1,...
            'delimiter','\t','emptyvalue', NaN,'headerlines',1);
% Reg ID string
    S=char(Reg_ID);

%% Loop through all the Reg_IDs looking for the profile you want
    for ii=1:length(S);index(ii)=strcmp(strtrim(S(ii,:)),profile2Analyse);end
        jj=find(index==1);
% Check that the profile exists in the file
    if isempty(jj)==1; % Profile missing in file
        disp(['Profile missing in survey recorded on ',datestr(dateNumber(mm))])
        qualityIndex(mm)=0;
        contourPosition(mm)=NaN;
    else % if exists then plot and extract contour
        OSGBx=easting(jj);
        OSGBy=northing(jj);
        xVar=chainage(jj);
        yVar=elevation(jj);
        [xVar,sortIndex]=sort(xVar);
        yVar=yVar(sortIndex);

        qualityIndex(mm)=1;
        data(nProfile).x=xVar;
        data(nProfile).z=yVar;
        data(nProfile).date=dateNumber(mm);
        nProfile=nProfile+1;
    end
        clear jj index

    end %% End file loop
        data(1).nProfiles=nProfile-1;


% put data on a regular grid
for ii=1:length(data)
         minx(ii)=min(data(ii).x);
         maxx(ii)=max(data(ii).x);
end
        xMin=min(minx);
        xMax=max(maxx);
        survey.x=xMin:xMax;
for ii=1:length(data)
        [data(ii).x,ia]=unique(data(ii).x);
        data(ii).z=data(ii).z(ia);
        survey.zArray(:,ii)=interp1(data(ii).x,data(ii).z,survey.x);
        survey.date(ii)=data(ii).date;
end

% Need to insure dates are in order
 [survey.date,iSort]=sort(survey.date);
 survey.zArray=survey.zArray(:,iSort);

 % Do timewise extrapolation
%  [n,m]=size(survey.zArray);
%  for ii=1:n
%     TimeLine=survey.zArray(ii,:);
%     TimeVal=survey.date;
%     clear jj
%     jj=find(isnan(TimeLine)==1);
%     TimeLine(jj)=[];
%     TimeVal(jj)=[];
%     if length(TimeLine)>4
%
%         survey.zArray(ii,:)=interp1(TimeVal,TimeLine,survey.date,'nearest','extrap');
%
%     end
%  end


   %% Plot Profiles
    fig1=figure;
    axes1 = axes('Parent',fig1);
    hold(axes1,'all');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);

    plot(survey.x,survey.zArray)
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    kk=find(qualityIndex==1);
    legend(datestr(dateNumber(kk)))
    grid on
    title(['Profile Reg ID = ' profile2Analyse])
    xlabel('Chainage (m)')
    ylabel('Beach profile height (m ODN)')
    h=gca;
    set(h,'FontSize',16);

%% Define verticle limits - users should adjust this to fit the area
 % zlevel defines the shoreline contour elevation (normally MHW)
 % Volume calculations are between zmin and zmax
    clc
    zmax=input('Please input the berm height (MSR/2 + 0.5 Hs,max) (m): ');
    zmin=input('Please input the mean water level (m): ');
    zlevel=(zmin+zmax) /2;
    % Auto-computez-level as the average of zmin and zmax
    survey.zlevel=zlevel;survey.zmax=zmax;survey.zmin=zmin;
     data(1).zLevel=zlevel;

      % Intergate wrt x to get volume
      for ii=1:length(data)%m
        xtemp=survey.x;
        ztemp=survey.zArray(:,ii);
        % Remove NaNs
        jNaN=find(isnan(ztemp));
        xtemp(jNaN)=[];
        ztemp(jNaN)=[];

        jj=find(ztemp>=zmax,1,'first');
        if isempty(jj);jj=1;end
        jj2=find(ztemp<=zmin,1,'first');
        if isempty(jj2);jj2=length(jj2);end
        index=jj:jj2;

        survey.volume(ii)=trapz(xtemp(index),ztemp(index)-min(ztemp(index)));
        survey.MCL(ii)=min(xtemp(index))+survey.volume(ii)/range(ztemp(index));
        xIndex=find(ztemp<=zlevel);xIndex=min(xIndex);
        xo(1)=xtemp(xIndex);xo(2)=xtemp(xIndex+1);
        yo(1)=ztemp(xIndex);yo(2)=ztemp(xIndex+1);
        [survey.shoreline(ii)]=myInterp(yo,xo,zlevel);
        clear xtemp ztemp
      end

%% Plot volume time-series
    fig2=figure;
    axes2 = axes('Parent',fig2);
    hold(axes2,'all');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    survey.z=nanmean(survey.zArray');
    survey.ProfileID=profile2Analyse;
    %
    yyaxis left;
    plot(survey.date,survey.volume,'k--','Marker','o'); hold on
    ax=gca;
    ax.YColor='black';
    ylabel('Volume [m^2]')
    %
    yyaxis right;
    plot(survey.date,survey.shoreline,'r--','Marker','o'); hold on
    plot(survey.date,survey.MCL,'g--','Marker','o'); datetick
%     ax=gca;
%     ax.YColor='black';
    ylabel('Shoreline Position [m] (red), MCL [m] (green)')
    xlabel('Year')
    grid on
    title(['Profile Reg ID = ' profile2Analyse])
    h=gca;
    set(h,'FontSize',16);

    %% Save output model shoreline
    disp('Save an ouput datafile - This file can be used to initialise ForCE or compare with ForCE output')
    [Filename,PathName,FilterIndex]=...
        uiputfile('*.mat','Input Filename for Output Model Shoreline',pwd);
    if FilterIndex==0
        disp('No data file selected - ending program.')
    else
        fileString=[PathName Filename];
        %surveyDate=survey.date;surveyVolume=survey.volume;
        save(fileString,'survey')
    end


return
%% End program

%% Function below...
function [y]=myInterp(xo,yo,x)

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
    y=(yo(1)*(xo(2)-x)+yo(2)*(x-xo(1)))/(xo(2)-xo(1));
end
