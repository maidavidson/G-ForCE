function [forecastDates, Var1Forecast,Var2Forecast,Var3Forecast]...
    =makingWavesFunction(dates,var1,var2,var3,startDate,trailYears,duration,dt)
%% Function or program name: makingWavesFunction
% Input
%
%   Inputs:
%       dates — input argument (description/units to be verified)
%       var1 — input argument (description/units to be verified)
%       var2 — input argument (description/units to be verified)
%       var3 — input argument (description/units to be verified)
%       startDate — input argument (description/units to be verified)
%       trailYears — input argument (description/units to be verified)
%       duration — input argument (description/units to be verified)
%       dt — input argument (description/units to be verified)
%
%   Outputs:
%       forecastDates — output argument (description/units to be verified)
%       Var1Forecast — output argument (description/units to be verified)
%       Var2Forecast — output argument (description/units to be verified)
%       Var3Forecast — output argument (description/units to be verified)
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

%% Input
% dates of the input wave data (matlab datenum format)
% var1 - var2 - var3: Input wave parameters, e.g. H, T, direction. NB make
% sure the input data is clean! No NaNs.
% startDate - matlab datenum value for simulation start date NB this should
% be within the timespan of the inmput wave data as the routine tags on
% measured anticeedent conditions before the forecast series.
% trailYears - number of trail years used for anticeededent condition
% calculation before the model start date.
% duration - duration of synthetic wave forecast  (days)
% dt - wave input sampling interval
%% Output
% ForecastDates - forecast dates
% Var1-3Forecast - forecasted parameters for input variables var1-3
% Written by Mark Davidson

%% Variables
 endDate=startDate+duration;
 forecastDates=startDate:dt:endDate;
 Var1Forecast=zeros(size(forecastDates));
 Var2Forecast=zeros(size(forecastDates));
 Var3Forecast=zeros(size(forecastDates));

%% Trail years
 trailTime=round(trailYears*365); % convert to round days
 iPatch=find(dates >= startDate-trailTime & dates < startDate);
 datesPatch=dates(iPatch);
 var1Patch=var1(iPatch);
 var2Patch=var2(iPatch);
 var3Patch=var3(iPatch);

%% Get rid of any possible months in the precition year (unbiased simulation)
 jj=find(dates >= startDate & dates <= endDate);
 dates(jj)=[];
 var1(jj)=[];
 var2(jj)=[];
 var3(jj)=[];

%% Get the year, month day etc of the wave measurements
 [Y,M,D,Hr,Min,Sec]=datevec(dates); Y=Y(:);

%% Strip out the different months from the wave pool and build an array for each month
 var1=var1(:);
 var2=var2(:);
 var3=var3(:);

for ii=1:12
   jj=find(M==ii);
   ArrayStr1(ii).month=['Month' num2str(ii) '_var1'];
   eval([char(ArrayStr1(ii).month) ' = [ii*ones(size(var1(jj))) Y(jj) var1(jj)];' ]);
   ArrayStr2(ii).month=['Month' num2str(ii) '_var2'];
   eval([char(ArrayStr2(ii).month) ' = [ii*ones(size(var2(jj))) Y(jj) var2(jj)];' ]);
   ArrayStr3(ii).month=['Month' num2str(ii) '_var3'];
   eval([char(ArrayStr3(ii).month) ' = [ii*ones(size(var3(jj))) Y(jj) var3(jj)];' ]);

% Randomise the month pools
    eval(['Years = ' char(ArrayStr1(ii).month) '(:,2);']);
    uniqueYears=unique(Years);
    for kk=1:length(uniqueYears)
        % find indexes of month pool data which have the current unique
        % year and month values.
        eval(['indexNeeded=find(' char(ArrayStr1(ii).month) '(:,2)==uniqueYears(kk));']);

       % Now asign the first column of the month pool with these index values a randon number
       randNo=rand;
       eval([char(ArrayStr1(ii).month) '(indexNeeded,1)=randNo;']);
    end
    % Strip out random numbers
    eval(['randColumn= ' char(ArrayStr1(ii).month) '(:,1);']);
    % now we can sort the month pool by the random index
    [dum,indexSort]=sort(randColumn);
    % Use this sorting index to sort the variables
    eval([char(ArrayStr1(ii).month) ' = '  char(ArrayStr1(ii).month) '(indexSort,:);']);
    eval([char(ArrayStr2(ii).month) ' = '  char(ArrayStr2(ii).month) '(indexSort,:);']);
    eval([char(ArrayStr3(ii).month) ' = '  char(ArrayStr3(ii).month) '(indexSort,:);']);
end

% Now create synthetic series
 [Y2,M2,D,Hr,Min,Sec]=datevec(forecastDates); Y2=Y2(:);M2=M2(:);
 YearMonth=[Y2 M2];
 UniqueYearMonth=unique(YearMonth,'rows');
 nMonths=length(UniqueYearMonth);
for ii=1:nMonths
    kk=find(YearMonth(:,1)==UniqueYearMonth(ii,1) & YearMonth(:,2)==UniqueYearMonth(ii,2));
    nPoints=length(kk);
    % Find a random starting point within the record
    monthNow=UniqueYearMonth(ii,2);
    eval(['nTot=length(' char(ArrayStr1(monthNow).month) ');'])
    startPoint=ceil((nTot-nPoints-1)*rand);
    endPoint=startPoint+nPoints-1;
    indexWanted=startPoint:endPoint;
    eval(['Var1Forecast(kk)=' char(ArrayStr1(monthNow).month) '(indexWanted,3);']);
    eval(['Var2Forecast(kk)=' char(ArrayStr2(monthNow).month) '(indexWanted,3);']);
    eval(['Var3Forecast(kk)=' char(ArrayStr3(monthNow).month) '(indexWanted,3);']);
end

%% Finally put on trailing years
forecastDates       = [datesPatch(:);forecastDates(:)];
Var1Forecast        = [var1Patch; Var1Forecast(:)];
Var2Forecast        = [var2Patch; Var2Forecast(:)];
Var3Forecast        = [var3Patch; Var3Forecast(:)];
