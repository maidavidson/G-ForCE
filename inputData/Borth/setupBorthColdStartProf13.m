% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
%% Setup file for Borth
% Mark Davidson June 2025

%% Input/output datafile names (base directory is the project name directory, in which two folders must exit: ontputData & inputData). Filenames can be user defined.
 filename.opFileName='./outputData/BorthOpFileCOLDStartProf14optFinalHunt.mat'; % output path/file name [NB o/p file of same name will be overwritten]
 filename.seaLevel='./inputData/Borth/seaLevelData0.mat'; % path/file containing sea level change data
 filename.waves='./inputData/Borth/BorthWavesShallowSynthetic.mat'; %
 filename.dissipationMap='./inputData/Borth/BorthDissipationMapProf13.mat'; %
 filename.tides='./inputData/Borth/BorthTidesHAwithlagCorr&Resid2.mat'; %
 filename.profileCS='./inputData/Borth/Borthprofile13LIDARspinUp.mat'; %
 filename.profileHS='./inputData/Borth/surveyData_09-Mar-2020_PROF13.mat'; %
 filename.NEB='NEB_Borth_P13_LIDAR.mat'; %
 filename.sourceSink=[];%'./inputData/Borth/sedimentSourceSinkBorthProfile13.mat'; %
 filename.executable='G_ForCE_v1_1.m';

%% Time-scales for model run, saving and plotting
model.startDate=[];%datenum(2020,01,01); % Model start date. if specified as [], then first profile date is used
model.NYears=1;%5.5;% Number of years of the model simulation
model.dateIncrement=1; % Saved data interval [days]. NB Bug known but unfixed for values<1
model.datePlotIncrement=7; % Plot interval

%% Key model setup parameters
model.dx=2; % Spatial resolution of modeled profile
model.dt=1/(6*24); % Model time-step in days
model.fact=0.5; % 1 - implicit, 0 - explicit, 0.5 Crank-Nicholson
model.spin=1; % 1 for spin-up (cold start), 0 for a hot start.
model.NSpin=1; % Number of years for model spin-up
model.avInterval=model.NSpin*365; % Number of days used to average dissipation map
model.MorphFB=1; % Morphodynamic Feedback to Hydro module. 0=none. 1=total.

%% Constants relating to the model environment
model.sedimentDensity=2650; % [kg/m^3]
model.waterDensity=1025; % [kg/m^3]
model.poreSpace=0.65; % [Dimensionless] (1-porosity)
model.D50(1)=0.048; % Sediment grain size [m]
model.D50(2)=0.24e-3; % Sediment grain size [m]
model.f=1; % factor involved in sediment threshold >=1 (default=1)

%% Key model performance parameters
model.k1(1)=  0.005;% Morphological response parameter - layer 1
model.k1(2)=  0.01;% Morphological response parameter - layer 2
model.k2   =  0.01; % Diffusion parameter

%% Tide correction parameters
model.tideCoeff=1.; % 1 = no correction
model.tideOffset=-0.0727; % 0 = no offset

%% Over-wash, through-flow and swash parameters
model.berm=4; % Berm height [m]
model.nOverwash=0; % Initialize over-wash parameter to zero
model.layers=2; % Number of sediment layers in the model
model.gravel=3; % Elevation (m ODN) of pure gravel
model.mixed=0.0;  % Elevation (m ODN) of mixed beach (landward limit of pure sand)
model.tideThreshold=inf; % height over which transport is +ve
model.IrThreshold=1.2;

%% Setup variables for wave dissipation model
model.alpha=1;
model.beta=0.1;
model.gamma=[];
model.ksw=0.025;

%% Swash model
model.swash=1;
model.swashFact=0.884;

%% Sea level rise algorithm parameters
model.transRes=2; % Sub-sampling resolution [m]
model.maxIteration=100; % Max no. of iterations
model.distSmoothed=5; % Region either side of gradient excess to be considered [m]
model.avelanche=0; % Initially set to zero. A value = 1 after model run indicates avelanching has occurred
model.avelancheOn=1;
model.maxGrad=0.35;
model.optMethod=0; % = 0 for zero net transport at both boundaries. = 1 for optimum sediment transport cost (e.g. sea wall)
model.BW=0; % Breakwater present. 0=No, 1=Yes
model.BWx=[]; % Cross-shore location of BW (must be specified for BW=1.

%% Constants
constant.days2seconds=24*60*60; % days to seconds conversion factor
constant.g=9.81; % [m/s^2] gravity
constant.k1=1/((model.sedimentDensity-model.waterDensity)*model.poreSpace*constant.g); % Sediment transport coefficient
