%% Function or program name: annimateForCEoutput
% User data
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
close all
clear
clc

%% User data
 disp('Select Directory for output file')
 selpath = uigetdir(path,'Select Dir for O/P file');
 opFile=[selpath,'/ModelRun2.avi'];
 FrameRate=10;
 dateStart=datenum(2000,1,1);
 dateFinish=datenum(2100,12,31);
 Res=30*8*1; % Extract every Res points [points]
 ymax=50;
 ymin=-15;
%   xmax=-1750;
%   xmin=-3000;

%% Load the model structure from the file
disp('LoadForCE ouput datafile')
 [file,path] = uigetfile('*.mat');
  if file~=0
   open(fullfile(path,file))
  end
 filename=fullfile(path,file);
 S=whos('-file',filename);
 SS=char(S.name);

%% load model and profile structure
load(fullfile(path,file),'model','profile','seaLevel')

%% load in profiles and dates for all years stored in the file
 kk=1;
 for ii =1:length(S)
   disp(char(S(ii).name));
   index(ii)=startsWith(SS(ii,:),'Profiles');
   if index(ii)==1
       variableName(kk,:)=SS(ii,:);
       mm=length(SS(1,:));
       dateVal=SS(ii,end-3:end);
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

%% Initialise elevation and time arrays
 zArray=[];
 tArray=[];

 %% Loop through all years extracting data into arrays
 for ii=jj(1):jj(end)
  variable1=['Profiles_',num2str(numDate(ii))];
  variable2=['Dates_',num2str(numDate(ii))];
  load(filename,variable1,variable2);
  eval(['zArray=[zArray,',variable1,'(:,1:Res:end)];']);
  eval(['tArray=[tArray,',variable2,'(1:Res:end)];']);
  eval(['clear ',variable1,' ',variable2])
 end

 %% Main program
 axis tight manual
 v=VideoWriter(opFile);
 v.FrameRate=FrameRate;
 open(v);
iCount=0;
jj=find(tArray>=dateStart & tArray<dateFinish);
m=length(jj);
zVal=  -profile.heq'+zArray;

% Find range of interest
 ix=find(-profile.heq>=ymin & -profile.heq<=ymax);
 xmin=profile.x(ix(1));
 xmax=profile.x(ix(end));

for ii=1:m
    iCount=iCount+1;
    SL=interp1(seaLevel.dates,seaLevel.z,tArray(jj(ii)),'linear','extrap');
    xSL=interp1(zVal(:,jj(ii)),profile.x,SL);
    xSL0=interp1(zVal(:,jj(ii)),profile.x,0);
    cla
    x=[profile.x,profile.x(end),profile.x(1)];
    y=[zVal(:,jj(ii));ymin;ymin];
    fill(x,y,'y')
    hold on
    plot(profile.x,zVal(:,1),'g:','LineWidth',3)
    plot([xmin xSL0],[0 0],'b:','LineWidth',2)
    plot([xmin xSL],[SL SL],'b','LineWidth',3)
    hold off
    xlim([xmin xmax])
    ylim([ymin ymax])
    legend('Modelled Profile','Original Profile','Original Mean Sea Level','Modelled Mean Sea Level',...
    'Location','northwest')
    grid on
    title(['Prediction date: ',datestr(tArray(jj(ii))),', Sea Level: ',num2str(SL),'m'])
    xlabel('Chainage[m]')
    ylabel('Elevation [m]')
    haxis=gca;
    set(haxis,'FontSize',16);
    drawnow

    frame=getframe(gcf);
    writeVideo(v,frame);

end

close(v)
