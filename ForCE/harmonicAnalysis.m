
%% Copyright notice
%   -----------------------------------------------------------------------
%   Copyright (C) 2014 Plymouth University
%
%       Mark Davidson Harmonic Analysis of Tides
%       Last revision: 20/11/2014
%
%       mdavidson@plymouth.ac.uk
%
%       Drake Circus
%       Plymouth, PL4 8AA
%
%       Input:  tidal elevation (eta)
%               Matlab datenum time (dnum)
%               flag=1 to plot results, otherwise  no plot
%       Output: Structure arays
%               HarmonicConstituents containing name, amplitude and phase
%               of components.
%               HarmonicData containing mean, rms error and Pierson coefficient
%               HarmonicPrediction containing time, predicted tide and
%               residuals
%
%   See <http://www.coastalprocesses.org/>.
%   ----------------------------------------------------------------------
%%
function [HarmonicConstituent,HarmonicPrediction,HarmonicData]=harmonicAnalysis(dnum,eta,flag)
%% Function or program name: harmonicAnalysis
% Copyright notice
%
%   Inputs:
%       dnum — input argument (description/units to be verified)
%       eta — input argument (description/units to be verified)
%       flag — input argument (description/units to be verified)
%
%   Outputs:
%       HarmonicConstituent — output argument (description/units to be verified)
%       HarmonicPrediction — output argument (description/units to be verified)
%       HarmonicData — output argument (description/units to be verified)
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

%HARMONICANALYSIS G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.

disp('Computing Harmonic analysis')

% if flag is not specified then don't plot data
if nargin==2;flag=0;end

C=eta; % The data array matrix C in course notes
dnum=dnum*24; % convert time in days to hours


%% Define 37 NOAA tidal radial frequencies (rf) (degrees/solar hour)
rf= ...
[28.9841042...
 30.0000000...
 28.4397297...
 15.0410690...
 57.9682083...
 13.9430351...
 86.9523163...
 44.0251732...
 60.0000000...
 57.4238319...
 28.5125828...
 90.0000000...
 27.9682083...
 27.8953552...
 16.1391010...
 29.4556255...
 15.0000000...
 14.4966936...
 15.5854435...
  0.5443747...
  0.0821373...
  0.0410686...
  1.0158958...
  1.0980331...
 13.4715147...
 13.3986607...
 29.9589329...
 30.0410671...
 12.8542862...
 14.9589310...
 31.0158958...
 43.4761581...
 29.5284786...
 42.9271393...
 30.0821381...
115.9364166...
 58.9841042];

% Define constituent names
Names=[...
'M2    ';...
'S2    ';...
'N2    ';...
'K1    ';...
'M4    ';...
'O1    ';...
'M6    ';...
'MK3   ';...
'S4    ';...
'MN4   ';...
'NU2   ';...
'S6    ';...
'MU2   ';...
'2N2   ';...
'OO1   ';...
'LAMDA2';...
'S1    ';...
'M1    ';...
'J1    ';...
'MM    ';...
'SSA   ';...
'SA    ';...
'MSF   ';...
'MF    ';...
'RHO1  ';...
'Q1    ';...
'T2    ';...
'R2    ';...
'2Q1   ';...
'P1    ';...
'2SM2  ';...
'M3    ';...
'L2    ';...
'2MK33 ';...
'K2    ';...
'M8    ';...
'MS4   '];

%% Evaluation of Harmonics...

% convert from degrees/hr to radians/hr
rf=(rf*2*pi)/360;

% Define number of coefficients and data points
M=length(C);
N=length(rf);

% Make A matrix (see lecture notes for definition of A)
A=ones(M,2*N+1);
jj=0;
for ii=1:2:2*N-1
    jj=jj+1;
    A(:,ii)  =sin(rf(jj)*dnum);
    A(:,ii+1)=cos(rf(jj)*dnum);
end

% Least squares to get unknown coefficients in matrix B
B=A\C;

% Reconstruct time series
Cm=A*B;

% Compute residual R
R=C-Cm;

% Define putput structure array
HarmonicPrediction.time=dnum./24;
HarmonicPrediction.tide=Cm;
HarmonicPrediction.residual=R;
HarmonicPrediction.measured=C;

% Errors and model skill
HarmonicData.RMSE=rms(R);
Pr=corrcoef(C,Cm);
HarmonicData.PiersonR=Pr(1,2);
HarmonicData.Mean=B(2*N+1);

%%
% Compute amplitude and phase and print to scrren
clc
fprintf('Harmonic Analysis Results \n')
fprintf('========================= \n')
fprintf('\n')
formatSpec='Root-mean-square error = %6.2f m, Pierson correlation coefficient = %4.2f \n';
fprintf(formatSpec,HarmonicData.RMSE,HarmonicData.PiersonR)
formatSpec='Mean tidal elevation   = %6.2f m \n';
fprintf(formatSpec,HarmonicData.Mean)
fprintf('\n')

formatSpec='Harmonic component: %4s, Period: %8.3f hrs, Amplitude: %6.3f metres, Phase: %8.3f degrees\n';
 jj=0;
for ii=1:2:2*N-1
    jj=jj+1;

    HarmonicConstituent(jj).Name        =Names(jj,:);
    HarmonicConstituent(jj).Period      =(2*pi)/rf(jj); % hrs
    HarmonicConstituent(jj).Amplitude   =sqrt(B(ii).^2 + B(ii+1).^2); % m
    HarmonicConstituent(jj).Phase       =atan(B(ii+1)/B(ii)); % radians


    fprintf(formatSpec,HarmonicConstituent(jj).Name,HarmonicConstituent(jj).Period,HarmonicConstituent(jj).Amplitude,(180/pi)*HarmonicConstituent(jj).Phase)
end

HarmonicData.B=B
HarmonicData.FormFactor = (HarmonicConstituent(4).Amplitude + HarmonicConstituent(6).Amplitude) ./ ...
                          (HarmonicConstituent(1).Amplitude + HarmonicConstituent(2).Amplitude);
if HarmonicData.FormFactor  <0.25; HarmonicData.Type='Semidiurnal';end
if HarmonicData.FormFactor >=0.25 & HarmonicData.FormFactor  <1.5; HarmonicData.Type='Mixed mainly semidiurnal';end
if HarmonicData.FormFactor >=1.5  & HarmonicData.FormFactor  <3.0; HarmonicData.Type='Mixed mainly diurnal';end
if HarmonicData.FormFactor  >3.0 ; HarmonicData.Type='Diurnal';end

fprintf('\n')
fprintf('Tidal Classification \n')
fprintf('==================== \n')
fprintf('\n')
formatSpec='Tidal Form Factor      = %6.2f \n';
fprintf(formatSpec,HarmonicData.FormFactor)
formatSpec='Tidal Classification   = %s  \n';
fprintf(formatSpec,HarmonicData.Type)


%%
% Plot comparison
if flag==1
subplot(2,1,1)
plot(dnum./24,C)
hold on
plot(dnum./24,Cm,'r')
title('Tidal displacement (m)')
xlabel('Date')
ylabel('metres')
legend('Measured','Predicted')
datetick('x','dd-mm-yyyy')
grid on
set(gca,'LineWidth',2,'FontSize',16)

% Residual plot
subplot(2,1,2)
plot(dnum./24,R)
datetick('x','dd-mm-yyyy')
title('Residual')
xlabel('Date')
ylabel('metres')
grid on
set(gca,'LineWidth',2,'FontSize',16)

end
disp('Done')
%%
