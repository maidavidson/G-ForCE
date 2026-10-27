function ftp_WWIII_forecast(~,~)
%% Function or program name: ftp_WWIII_forecast
% FTP_WWIII_FORECAST G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       ~ — input argument (description/units to be verified)
%       ~ — input argument (description/units to be verified)
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

%FTP_WWIII_FORECAST G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.

clc
clear all
close all
disp(['Last run at: ', datestr(now)])

%% Generate directories to be searched
for ii=1:8
    searchDir(ii)=cellstr(['multi_1.',datestr(today-8+ii,'yyyymmdd')]);
    folderDate(ii)=cellstr([datestr(today-8+ii,'yyyymmdd')]);
    folderDnum(ii)=datenum(char(folderDate(ii)),'yyyymmdd');
end

if exist('waveFile.mat') ==2
    disp('Existing wave file found and loaded')
    load waveFile
else
    % define wave parameters here
    % Define o/p array
    wave.N=1;
    wave.H=[];
    wave.T=[];
    wave.Dir=[];
    wave.date=[];
    wave.flag={};
    wind.Speed=[];
    wind.Dir=[];
    wave.LastFolderDate=folderDnum(1);
end


wave.plot =1; % 0 = no plot

%% delete all forecast record
jj = find(~cellfun('isempty',strfind(wave.flag,'forecast')));jj=jj';
wave.H(jj,:)=[];
wave.T(jj,:)=[];
wave.Dir(jj,:)=[];
wave.flag(jj)=[];
wave.date(jj)=[];
wave.N=length(wave.H)+1;
wind.Speed(jj,:)=[];
wind.Dir(jj,:)=[];

%% Now see what drectories need to be searched...
    dateLastFile=datenum(char(wave.LastFolderDate),'yyyymmdd');
    jj=find(folderDnum<dateLastFile);
    searchDir(jj)=[];
    folderDate(jj)=[];
    folderDnum(jj)=[];

%% Locations
% Sevenstones
latSS=50.102;lonSS=360-6.100;
latPP=50.35276;lonPP=360-5.17503;
latPL=50.06223;lonPL=360-5.30801;
wave.location=[latSS,lonSS;latPP,lonPP;latPL,lonPL];
wave.site={'Sevenstones','Perranporth','Porthlevan'};

%% Download the last 7 days data
for ii=1:length(searchDir)

%% Open ftp site
ftpobj = ftp('ftpprd.ncep.noaa.gov');

% Change dir to archive directory
baseFolder=['pub/data/nccf/com/wave/prod/',char(searchDir(ii))];
NewFolder=cd(ftpobj,baseFolder);

% Retrieve the file names and count them
files=dir(ftpobj,'*multi_1.glo_30m.t*f000*.grib2');
Nfiles=length(files);

%% Loop through all files now downloading gribs and extracting data
for jj=1:Nfiles
    mget(ftpobj,char(files(jj).name));
    wave.flag(wave.N)={'nowcast'};
    disp(['Extracting data from: ',char(files(jj).name) ])
    [wave,wind]=extractGribData(char(folderDate(ii)),char(files(jj).name),wave,wind);
    LastFile=char(files(jj).name);

end % nowcast files in each dir

%% If on last day (ii==8), then get the forecast
if ii==length(searchDir)
    files=dir(ftpobj,[LastFile(1:22),'*.grib2']);
    Nfiles=length(files);
    for jj=1:3:Nfiles
        try
        mget(ftpobj,char(files(jj).name));
        wave.flag(wave.N)={'forecast'};
        char(files(jj).name);
        disp(['Extracting data from: ',char(files(jj).name) ])
        [wave,wind]=extractGribData(char(folderDate(ii)),char(files(jj).name),wave,wind);
        catch
           disp('Could not open a connection to "ftpprd.ncep.noaa.gov", port "21".')
           wave.N=wave.N-1;
        end
  end % forecast files in each dir
end   % if

% close ftp connection and turn command line echo off
close(ftpobj)
delete *.grb2 *.gbx9 *.grib2

end % end directories


%% Clean up and save data
    % delete non-unique records
    [wave.date,ia,ic]=unique(wave.date);
    wave.H=wave.H(ia,:);
    wave.T=wave.T(ia,:);
    wave.Dir=wave.Dir(ia,:);
    wave.flag=wave.flag(ia);
    wave.N=length(wave.H)+1;
    wind.Speed=wind.Speed(ia,:);
    wind.Dir=wind.Dir(ia,:);

    delete *.grb2 *.gbx9 *.grib2 *.ncx
    delete(['/Users/markdavidson/Documents/MATLAB/*.ncx'])
    wave.LastFile=LastFile;
    wave.LastFolderDate=folderDate(end);
    cd wave_data
    save waveFile wave wind
    cd ../
%% Plot wave data
plotWaves(wave,wind)

%% End programmme
%%
function plotWaves(wave,wind)
%PLOTWAVES G_ForCE support routine; see function signature for inputs and outputs.
    fig1=figure;
    axes1 = axes('Parent',fig1);
    hold(axes1,'all');

    %% Maximise
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    %% Plot wave height
    subplot(5,1,1)
    plot(wave.date,wave.H)
    hold on
    plot([now, now],[0,max(6.64+0.5,max(max(wave.H)))],'k--')
    plot([wave.date(1),wave.date(end)],[6.64, 6.64],'r-.')
    title(['Significant Wave Height last updated: ',datestr(now)])
    ylabel('[m]')
    xlabel('Date')
    xticks(floor(today-7): ceil(today+8))
    datetick('x','dd/mm','keepticks')
    %xlim(floor(today-7): ceil(today+7))
    grid on
    legend(wave.site,'Location','northwest')
    h=gca;
    set(h,'FontSize',14);
    %% plot wave period
    subplot(5,1,2)
    plot(wave.date,wave.T)
    hold on
    plot([now, now],[0,max(max(wave.T))],'k--')
    title('Peak Wave Period')
    ylabel('[s]')
    xlabel('Date')
    xticks(floor(today-7): ceil(today+8))
    datetick('x','dd/mm','keepticks')
    %xlim(floor(today-7): ceil(today+7))
    grid on
    h=gca;
    set(h,'FontSize',14);
    %% Plot wave direction
    subplot(5,1,3)
    plot(wave.date,wave.Dir)
    hold on
    plot([now, now],[0,360],'k--')
    plot([wave.date(1),wave.date(end)],[180, 180],'k-.')
    plot([wave.date(1),wave.date(end)],[90, 90],'k-.')
    plot([wave.date(1),wave.date(end)],[270, 270],'k-.')
    title('Peak Wave Direction')
    ylabel('[Degrees]')
    xlabel('Date')
    xticks(floor(today-7): ceil(today+8))
    datetick('x','dd/mm','keepticks')
    %xlim([floor(today-7) ceil(today+7)])
    ylim([0,360])
    grid on
    h=gca;
    set(h,'FontSize',14);
    %% Plot wind speed
    subplot(5,1,4)
    plot(wave.date,wind.Speed)
    hold on
    plot([now, now],[0,360],'k--')
    hold on
    plot([wave.date(1),wave.date(end)],[17, 17],'r-.')
    plot([wave.date(1),wave.date(end)],[8, 8],'g-.')
    title('Wind Speed')
    ylabel('[m/s]')
    xlabel('Date')
    xticks(floor(today-7): ceil(today+8))
    datetick('x','dd/mm','keepticks')
    %xlim([floor(today-7), ceil(today+7)])
    ylim([0 max(max(wind.Speed))])
    grid on
    h=gca;
    set(h,'FontSize',14);
    %% Plot wind direction
    subplot(5,1,5)
    plot(wave.date,wind.Dir)
    hold on
    plot([now, now],[0,360],'k--')
    plot([wave.date(1),wave.date(end)],[180, 180],'k-.')
    plot([wave.date(1),wave.date(end)],[90, 90],'k-.')
    plot([wave.date(1),wave.date(end)],[270, 270],'k-.')
    title('Wind Direction')
    ylabel('[Degrees]')
    xlabel('Date')
    %xlim([floor(today-7), ceil(today+7)])
    xticks(floor(today-7): ceil(today+8))
    datetick('x','dd/mm','keepticks')
    ylim([0,360])
    h=gca;
    set(h,'FontSize',14);

    %% save fig as image
    saveas(fig1,'./output/LatestWave&WindForecast.png')

disp(['All done at: ', datestr(now)])
end

%%
function [wave,wind]=extractGribData(folderDate,fileName,wave,wind);
%EXTRACTGRIBDATA G_ForCE support routine; see function signature for inputs and outputs.

%% Now do plotting and data extraction
    nco=ncgeodataset(fileName);
    echo off all


% Extract the Significant Height of Swell Waves field.
    param={ 'Significant_height_of_combined_wind_waves_and_swell_surface',...
            'Primary_wave_mean_period_surface'...
            'Primary_wave_direction_degree_true_surface',...
            'Wind_direction_from_which_blowing_degree_true_surface',...
            'Wind_speed_surface'};
    lat=nco{'lat'}(:);
    lon=nco{'lon'}(:);

% Forecast horizon
    tVal=str2double(fileName(18:19)) + str2double(fileName(23:25));
    n=length(tVal);

   % Work out date and time
    dnum(wave.N)=datenum(folderDate,'yyyymmdd') + tVal/24;


    waveheight=nco{param(1)}(1,:,:);
    waveheight=double(squeeze(waveheight));
    lat=double(lat);
    lon=double(lon);
    waveperiod=nco{param(2)}(1,:,:);
    waveperiod=double(squeeze(waveperiod));
    wavedirection=nco{param(3)}(1,:,:);
    wavedirection=double(squeeze(wavedirection));
    windspeed=nco{param(5)}(1,:,:);
    windspeed=double(squeeze(windspeed));
    winddirection=nco{param(4)}(1,:,:);
    winddirection=double(squeeze(winddirection));

    if wave.plot ==1
        m_proj('miller','lat',[min(lat(:)) max(lat(:))],...
        'lon',[min(lon(:)) max(lon(:))])
        m_pcolor(lon,lat,waveheight);
        shading flat;
        m_coast('patch',[.7 .7 .7]);
        m_grid('box','fancy')
        colorbar
        title(['WAVEWATCH III Significant Wave Height. Forecast: ', datestr(double(dnum(wave.N)))]);
        drawnow
    end
    % Record variables
    [nn,m]=size(wave.location);
    for kk=1:nn
        wave.date(wave.N)=double(dnum(wave.N));
        wave.H(wave.N,kk)=interp2(lon,lat,waveheight,wave.location(kk,2),wave.location(kk,1));
        wave.T(wave.N,kk)=interp2(lon,lat,waveperiod,wave.location(kk,2),wave.location(kk,1));
        wave.Dir(wave.N,kk)=interp2(lon,lat,wavedirection,wave.location(kk,2),wave.location(kk,1));
        wind.Speed(wave.N,kk)=interp2(lon,lat,windspeed,wave.location(kk,2),wave.location(kk,1));
        wind.Dir(wave.N,kk)=interp2(lon,lat,winddirection,wave.location(kk,2),wave.location(kk,1));
    end
    wave.N=wave.N+1;

end % subroutine

end
