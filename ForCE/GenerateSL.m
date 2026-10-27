%% Function or program name: GenerateSL
% Predict sea level rise curve
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

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% Predict sea level rise curve
clear
close all
clc

% Formula = z = a + bt + c t^2

% known points
 a=0; % Offset
 b=0.002 / 365.25; %m/day - linear trend
 startYear=input('Please input the year when sea level started to rise, e.g. 2025.  ');
 knowndZ=input('Please input the predicted sea level rise 100 years from the start date, e.g. 0.5 m.  '); % Known dz
 knownTime=100; % Known time for known dz

 % Convert time to days
 years2days=365.25;
 knownTime=knownTime*years2days;

% Compute c
c= (1/(knownTime)) * ( ((knowndZ-a)/knownTime)  - b); % [m/day^2] Acceletaion term

% Compute timebase for calcs
startDate=datenum(1980,1,1);
startRise=datenum(startYear,1,1);
dt=1; % Sampling rate [days]
nYears=320;
nDays=nYears*365.25;
endDate=startDate+nDays;
seaLevel.dates=startDate:dt:endDate;

% Compute sea level rise
jj=find(seaLevel.dates<startRise);
N=length(seaLevel.dates);
kk=jj(end):N;
seaLevel.z(jj)=a;
seaLevel.z(kk)=a + b* (seaLevel.dates(kk)-startRise) + c *(seaLevel.dates(kk)-startRise).^2 ;

% Plot
plot(seaLevel.dates,seaLevel.z,'LineWidth',3)
datetick('x','yyyy')
ylabel('Sea level Rise [m]')
xlabel('Year')
grid on
box on
set(gca,'LineWidth',3,'FontSize',16)

% Save Data
saveFile=input('Do you wish to save the file (y/n)? ','s');
if saveFile=='y' | saveFile=='Y'
    % get input filename
    [FileName,PathName] = uiputfile({'*.mat'},'Save sea level file as...');
    % Inputs data
    outputFile=[PathName FileName];
    save(outputFile, 'seaLevel');
    disp('File saved')
else
    disp('File not saved')
end

disp('Done')
