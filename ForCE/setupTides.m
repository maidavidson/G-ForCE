%% Function or program name: setupTides
% Harmonic analysis of measured or modeled data to fill gaps and extend
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
% Harmonic analysis of measured or modeled data to fill gaps and extend
% tidal records - Mark Davidson 2025


close all
clear all

% Synthetic wave generator (Mark Davidson 2025)
disp('Harmonic tide generator')

% Get waves file
 disp('Please load the example tide file containg modeled or measured waves')
 disp('This file should be a .mat format containing the structure tide with the fields dates and z (m AOD)')
 [file,path] = uigetfile('*.mat');

% Load data
 load(fullfile(path,file))

 % define period of simulation
        dnumStart   =datenum(1980,1,1);
        dnumEnd     =now+(100*365);
        dnum        =24*(dnumStart:tide.dt:dnumEnd)';

% Pass the data in the structure tide for harmonic analysis
flag=1; % flag equals 1 for plotting the results
[HarmonicConstituent,HarmonicPrediction,HarmonicData]=harmonicAnalysis(tide.dates,tide.z,flag);

% Define B-Array coefficients
B=HarmonicData.B;
[MM,NN]=size(HarmonicConstituent);
for ii=1:NN; T(ii)=HarmonicConstituent(ii).Period; end
rf=2*pi./T;

% Define number of coefficients and data points
M=length(dnum);
N=length(rf);

% Make A matrix (see lecture notes for definition of A)
A=ones(M,2*N+1);
jj=0;
for ii=1:2:2*N-1
    jj=jj+1;
    A(:,ii)  =sin(rf(jj)*dnum);
    A(:,ii+1)=cos(rf(jj)*dnum);
end

C=A*B;

tideOld=tide;
tide.z=C;
tide.dates=dnum./24;
% Compute residuals and remove outliers
tide.residual=zeros(size(tide.dates));
jj=find(tide.dates >= tideOld.dates(1) & tide.dates <= tideOld.dates(end));
z=interp1(tideOld.dates,tideOld.z, tide.dates(jj));
tideRes=z-tide.z(jj);
[~,ir]=rmoutliers(tideRes);
tideRes(ir)=0;
tide.residual(jj)=tideRes;

%%
figure
plot(tide.dates,tide.z+tide.residual)
hold on
plot(tide.dates,tide.z)
title('Tidal displacement (m)')
xlabel('Date')
ylabel('metres')
legend('Astronomic tide + surge','Astronomic tide')
datetick('x','dd-mm-yyyy')
grid on
set(gca,'LineWidth',2,'FontSize',16)


% Save Data
saveFile=input('Do you wish to save the file (y/n)? ','s');

if saveFile=='y' | saveFile=='Y'
    % get input filename
    [FileName,PathName] = uiputfile({'*.mat'},'Save as (NB Give the station location in the Name');
    % Inputs data
    outputFile=[PathName FileName];
    matlabTime=dnum/24;tideHeight=C;
    save(outputFile, 'tide');
    disp('File saved')
else
    disp('File not saved')
end

disp('Done')
