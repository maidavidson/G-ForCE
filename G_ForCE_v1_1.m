%% Function or program name: G_ForCE_v1_1
% G-ForCE
%
%   Inputs:
%       Workspace file reference: filename.seaLevel
%       Workspace file reference: filename.waves
%       Workspace file reference: filename.profileCS
%       Workspace file reference: filename.tides
%       Workspace file reference: filename.dissipationMap
%       Workspace file reference: filename.sourceSink
%       Workspace file reference: filename.NEB
%       Workspace file reference: filename.profileHS
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
% G-ForCE
% Gravel - Forecasting Coastal Evolution model.
% Based on the ForCE model for sandy coasts.
% Innovations:
%   - Optional two interactive layers of sand and gravel, plus a mixed zone
%   - improved wave model including the swash zone
%   - overwash and barrier through flow included
% Written by Mark Davidson
% NERC Gravel Beach Project
% 2025-2026


% Program started
 disp(['Program started at: ',string(datetime('now'))])
 tic

% Setup paths. Modify first one with GUI later
 addpath(path)
 addpath('./ForCE/')
 addpath('./linear/')

% Load data
load(filename.seaLevel)
load(filename.waves)
load(filename.profileCS)
load(filename.tides)

% Compute wave Power from waves file.
disp('Computing wave power')
[k]=kLinear(wave.Tp,wave.depth);
[wave.n]=nLinear(k,wave.depth);
[wave.c]=cLinear(wave.Tp,wave.depth);
[wave.E]=ELinear(wave.H);
[wave.P]=PLinear(wave.E,wave.c,wave.n);
 wave.Pmean=mean(wave.P);
disp('Finished computing wave power')

% Change wave angle from bearing to rel. to normal
if isfield(model,'meanBearing')
 [wave.Meanbearing, wave.theta]=convertAngle(wave.Dir,model.meanBearing);
else
 [wave.Meanbearing, wave.theta]=convertAngle(wave.Dir);
end

% Save empty output file to be appended to later
 save(filename.opFileName);

% compute DOC
 [model.DOC,~,~]=Hallermeier(wave.H,wave.Tp,wave.dt);
  model.DOC=model.DOC-min(tide.z); % Add in tide

% Add residual to tide data and apply corrections if required
if isfield(model,'tideOffset')
 tide.z=model.tideCoeff*(tide.z+tide.residual)+ model.tideOffset;
else
  tide.z=tide.z+tide.residual;
end

% Load profile and put in ForCE format with x inc. positively onshore.
% Extrapolate to DOC
if model.spin==0 && exist(filename.dissipationMap,'file')==2 % If a hot start then re-interpolate starting profile onto the original spin-up chainage
  load(filename.dissipationMap,'-mat')
  htemp=interp1(-flip(survey.x),-flip(survey.z),x);
  jj=find(~isnan(htemp));
  h0=heq;
  h0(jj)=htemp(jj);h0=smooth(h0,3);h0=h0(:)';
  clear jj htemp

else % If model spinning up then extrapolate the spin-up profile
    [x,zout]=extrapolateSurvey(survey.x,survey.z,-(model.DOC+1),0,model.dx); % extrapolates short profiles to DOC
    x=x(:)';zout=zout';
    h0=-flip(zout);h0=smooth(h0,3); h0=h0(:)'; % h0 = initial profile at t=0
    x  =-flip(x);
end

% Record model x-values
 model.x=x;

% Use the starting profile date to start the simulation
if isempty(model.startDate)
   model.startDate=survey.date(1);
end


% Setup initial sea level change parameters
 dxSLRsum=0; dxSLR=0; SLROld=0;
 SLR0=interp1(seaLevel.dates, seaLevel.z,model.startDate,'linear','extrap');

% Load sediment source-sink file if it exists?
if exist(filename.sourceSink,'file') ==2 % file already exists
  load(filename.sourceSink)
  disp('Source/sink file loaded')
  Vold=interp1(sourceSink.dates, sourceSink.V,model.startDate,'linear','extrap');
else
   Vold=0;
   sourceSink.V=[0 0];
   sourceSink.dates=[wave.dates(1) wave.dates(end)];
   sourceSink.units='m^3/m coast';
end

% Load non-erodible bed file if it exists?
if exist(filename.NEB,'file')
    load(filename.NEB)
    xNEB=NEB.x(:)';zNEB=NEB.z(:)';
    hNEB=-flip(zNEB); hNEB=hNEB(:)';
    xNEB=-flip(xNEB);
    hNEB=interp1(xNEB,hNEB,x,'linear','extrap');
    xNEB=x;
else
    NEB.x=[min(survey.x) max(survey.x)];
    NEB.z=[-40 -40];
    NEB.date=model.startDate;
    NEB.SW=0;
    NEB.xSW=[];
end

% Compute beach gradient
  model.beachGradient=abs(gradient(-h0,x));

% Compute estimate of Berm height (used in SL translation model)
if ~isfield(model,'berm')
    [bermHeight]=bermCalculator(wave,tide,model,SLR0);
else
    bermHeight=model.berm;
end

% Test to see if spin-up is required?
if model.spin==1
    nreps=2;
else
    nreps=1;
end

%% Cold/Hot-start counter/loop%%%%%%%%%%%%%%%%%
for irep=1:nreps

% Compute initial mean dissipation (spin-up) or load file (no-spin-up)
if model.spin~=1 && exist(filename.dissipationMap,'file')==2
    disp('Hot-start')
    load(filename.dissipationMap)
    % Storm cut algorithm here to modify equilibrium profile at run-up
    % limit
    heq0=heq; % Initial equilibrium shoreline
    if isfield(model,'stormCut')  && model.k2(1)~=0;[heq]=stormCut(x,heq,tide,model);end
    % Load in the hot-start profile now
    load(filename.profileHS)
    htemp=interp1(-flip(survey.x),-flip(survey.z),x);
    jj=find(~isnan(htemp));
    h0=heq;
    h0(jj)=htemp(jj);h0=smooth(h0,3);h0=h0(:)';
    clear jj htemp

    % Apply NEB to h0 so that it isn't passed to zPrime
    h0=min(h0,hNEB);
    % Setup zPrime so that the initial profile is given by h0=zPrime-heq
    zPrime=heq-h0;
    % re-set spin-up indicator to hot start
    model.spin=0;
else
    disp('Model spin-up in process')
    disp('Computing initial dissipation map')
    heq=h0;heq0=heq;
    [Deq,Deq2,depthAttenuation,overwash,maxRunup_x,maxWaterDepth,Gradient,swashIndex,Ir]=findMeanDissipation_3(x,heq,model,wave,tide,seaLevel,model.startDate,model.avInterval);
    if isfield(model,'stormCut') && model.k2(1)~=0;[heq]=stormCut(x,heq,tide,model);end
    model.meanGradient=Gradient;
    save(filename.dissipationMap,"Deq","Deq2","heq","x")
    disp('Finished computing initial dissipation map')
    zPrime=heq-h0;
end

% If run in SLR mode
if model.k1(1)==0
    heq0=h0; heq=h0; zPrime=0*zPrime;
end

% Record profile information
profile.zeq=-heq;
profile.zeq0=-heq0;
profile.z0=-h0;
profile.h0=h0;
profile.heq=heq;
profile.x=x;
if isfield(survey,'RegID')
 profile.ID=survey.RegID;
else
 profile.ID=survey.ProfileID;
end
profile.t0=model.startDate;

% Update model structure with equilibrium dissipation map
model.Deq=Deq;
model.Deq2=Deq2;
model.x=x;

% Now do the profile evolution. Start by defining model dates
if model.spin==0
    mdates=model.startDate:model.dt:model.startDate+365*model.NYears;
    model.Npts=length(mdates);
else
    mdates=model.startDate:model.dt:model.startDate+model.avInterval;
    model.Npts=length(mdates);
    zPrimeArray=zeros(length(x),model.Npts);
    DissArray=zeros(length(x),model.Npts);
end

% Define triggers and counters
jj=0; dateStore=0;
jRock=[];

% Save data so far [add model structure with model variables]
 save(filename.opFileName,'wave', 'profile', 'tide','seaLevel','model', '-append')


% Setup plot interval
 % An infinite plotting interval disables time-step plots entirely.
% This does not affect the independent data-saving interval.
if isinf(model.datePlotIncrement)
    nextPlot = inf;
else
    nextPlot = model.startDate;
end

% Initialise max run-up
 model.maxRunupMax_x(1)=-inf;
 model.maxRunupMax_z(1)=-inf;

% No. points in fortnight
 Nday=round(1/model.dt);

% Show lightweight Command Window progress only when time-step plots are off.
showPercentProgress = isinf(model.datePlotIncrement);
lastPercent = -1;
if showPercentProgress
    if model.spin == 1
        progressPhase = 'Cold start';
    else
        progressPhase = 'Hot start';
    end
    fprintf('%s: %3d%% complete', progressPhase, 0);
    lastPercent = 0;
end

%%%% BEGIN MAIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:model.Npts

% Record date
  model.dates(i)=mdates(i);

% % Re-evaluate recent berm height
 if i>Nday
     model.berm(i)=max(model.maxRunup_z(end-Nday+1:end));
     model.berm(i)=min(model.berm(i), max(-h0'));
 else
     model.berm(i)=bermHeight;
 end

% Do sea level adjustment
 heqOld=heq; % record heq prior to adjustment
 zPrimeOld=zPrime;
 dxSLROld=dxSLR;


 [heq,Vold,model,dxSLR,SLR,SLROld,flagUpdate]=SLR_Function(model,seaLevel,sourceSink,x,heq0,heq,mdates(i),Vold,dxSLR,SLR0,SLROld);
 if flagUpdate==1 % If sea level adjustment is made then do this...
      disp('SL/flux update complete')
      if isfield(model,'stormCut') && model.k2(1)~=0;[heq]=stormCut(x,heq,tide,model);end
      dxSLRsum=dxSLR; % Total profile displacement

      %if model.k2(1)~=0; zPrime=-heqOld+heq+zPrimeOld; clear zPrimeOld heqOld;end % Update zPrime only if diffusion is activated

      % Translate perturbation inline with equilibrium profile shift
      zPrime=interp1(x+(dxSLR-dxSLROld),zPrimeOld,x,'Linear','extrap');

      % avelanche if required
      maxGradOld=model.maxGrad;
      model.maxGrad=max(abs(gradient(-heq0,x)));
      [zout,~]=avelanche(x,-heq,model);
      model.maxGrad=maxGradOld; clear maxGradOld
      heq=-zout;
 end

 % Implement non-erodible bed
 heq=min(heq,hNEB);

% Call the ForCE algorithm
if model.k1(1)~=0 || model.k2(1)~=0

 [zPrime,q,D,DSand,model,maxRunup_x,maxWaterDepth,model.tanb(i),swashIndex,model.Ir(i),model.overwash(i)]=...
     ForCEfunction3(x,heq,zPrime,model,constant,wave,tide,seaLevel,mdates(i),jRock,dxSLRsum,SLR,nreps);

    if isempty(swashIndex)
        model.swashIndexMin(i)= NaN;
        model.swashIndexMax(i)= NaN;
        model.maxRunup_x(i)=maxRunup_x;
    else
        model.swashIndexMin(i)= min(swashIndex);
        model.swashIndexMax(i)= max(swashIndex);
        model.maxRunup_x(i)=x(swashIndex(end));
    end

    if i>1; model.maxRunupMax_x(i)=max(maxRunup_x,model.maxRunupMax_x(i-1));model.maxRunupMax_z(i)=max(maxWaterDepth,model.maxRunupMax_z(i-1));model.maxRunupMax_z(i)=min(model.maxRunupMax_z(i),max(-heq+zPrime) ) ;end
    model.maxRunup_z(i)=maxWaterDepth;
end


% avelanche if required
  [zout,iterations]=avelanche(x,-heq+zPrime,model);
  zPrime=zout+heq;

if iterations >0
  disp('avelanche')
end

%% Program output here
% Assign new output variable names based on year
 [thisYear,MM,DD,hh,mm,ss]=datevec(mdates(i));
 clear MM DD hh mm ss
 if i==1
  nextYear=thisYear+1;
  variable1=['Profiles_',num2str(thisYear)];
  variable2=['Dates_',num2str(thisYear)];
  variable3=['Diss_',num2str(thisYear)];
  iFile=1;
 end

 if thisYear==nextYear
  nextYear=nextYear+1;
  save(filename.opFileName,variable1, variable2,variable3, '-append')
  eval(['clear ',variable1,' ', variable2])
  variable1=['Profiles_',num2str(thisYear)];
  variable2=['Dates_',num2str(thisYear)];
  variable3=['Diss_',num2str(thisYear)];
  iFile=1;
 else
  variable1=['Profiles_',num2str(thisYear)];
  variable2=['Dates_',num2str(thisYear)];
  variable3=['Diss_',num2str(thisYear)];
 end

if model.spin==1
    zPrimeOut= -(heq-zPrime);
    zPrimeArray(:,i)=zPrimeOut(:);
end

% Assign updated values to arrays for output
% Update model output and plot
if dateStore<=mdates(i)
   dateStore=mdates(i)+model.dateIncrement;
   % Need to incorporate changing equilibrium profile here and add
   % this to hPrime output
   zPrimeOut= -(heq-zPrime);
if exist(filename.NEB,'file')
   jRock=find(zPrimeOut(:)<-hNEB(:));
   zPrimeOut(jRock)=-hNEB(jRock);
end
   eval([variable1,'(:,iFile)= zPrimeOut(:);']);
   eval([variable2,'(iFile)= mdates(i);']);
   eval([variable3,'(:,iFile)= D(:);']);
   iFile=iFile+1;
end % End if time to store data

%% plot results if time to do so
if mdates(i)>=nextPlot
   nextPlot=nextPlot+ model.datePlotIncrement;
clf
plot(x,-heq,'g--','LineWidth',2);
hold on
plot(x,zPrimeOut,'k','LineWidth',3);
plot(x,-h0,'r','LineWidth',1);
if exist(filename.NEB,'file')
    plot(xNEB,-hNEB,'b.','LineWidth',4);
end
plot(x,zPrime,'r','LineWidth',3)
if model.k1(1)~=0;yline(model.maxRunupMax_z(i),'b','LineWidth',3,'LineStyle',':');end
%plot(model.maxRunupMax_x(i),model.maxRunupMax_z(i),'r*','MarkerSize',12)
ylim([min(tide.z),max(1.1*max(-h0),max(-heq+zPrime))])
jx=find(zPrimeOut>=-2,1,'first');
xlim([x(jx), x(end)])
grid on
box on
set(gca,'LineWidth',3,'FontSize',24)
title([datestr(model.dates(i)),'. Sea level rise = ',num2str(round(1000*SLR)/1000),' [m]. Profile adjustment = ',num2str(round(100*dxSLRsum)/100),' [m]' ])
xlabel('Chainage [m]')
ylabel('Elevation [m AOD])')
drawnow
end


% Report progress at whole-percentage increments (at most 101 updates).
if showPercentProgress
    percentComplete = floor(100 * i / model.Npts);
    if percentComplete > lastPercent
        fprintf('\r%s: %3d%% complete', progressPhase, percentComplete);
        lastPercent = percentComplete;
    end
end

end %%%%%%%%%%%%%%%%%%% END MAIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if showPercentProgress
    fprintf('\n');
end


%% Compute equilibrium conditions in spin-up cycle
if model.spin==1
    heq=-mean(zPrimeArray,2);heq=heq';
    disp('Computing final dissipation map')
    [Deq,Deq2,depthAttenuation,overwash,maxRunup_x,maxWaterDepth,Gradient,swashIndex,Ir]=findMeanDissipation_3(x,heq,model,wave,tide,seaLevel,model.startDate,model.avInterval);
    model.meanGradient=Gradient;
    save (filename.dissipationMap, "Deq","Deq2" ,"heq","x")
    disp('Finished computing final dissipation map and equilibrium profile')
    model.spin=0;
    model.nOverwash=0;
end

% Save the final year's accumulated output after its last time step.
% This does not change the numerical model or the annual rollover logic.
if model.spin==0
    save(filename.opFileName,variable1,variable2,variable3,'-append')
end

end %%%%%%%%%%%%%%%%%%% END SPIN-UP & REPEAT LOOP %%%%%%%%%%%%%%%%%%%%%%%%%

% Save final values
 save(filename.opFileName,'wave', 'profile', 'tide','seaLevel','model', '-append')

%% Now load output data and plot
fname=filename.opFileName;
dateStart=datenum(1980,1,1);
dateFinish=datenum(2100,1,1);
[tArray,zArray]=getModelOutput(fname,dateStart,dateFinish);

figure
subplot(2,1,1)
plot(x,zArray,'k','LineWidth',3);
hold on
plot(x,-heq0,'w.','LineWidth',3);
if exist(filename.NEB,'file')
 plot(xNEB,-hNEB,'r.','LineWidth',6);
end
title(['Profile ID: ',profile.ID], 'Interpreter', 'none')
ylim([-model.DOC,1.1*max(-h0)])
jx=find(heq<=model.DOC,1,'first');
xlim([x(jx) x(end)])
grid on
box on
set(gca,'LineWidth',3,'FontSize',24)

subplot(2,1,2)
zplot=zArray'+heq; % plot change rel. to starting profile
jj=find(abs(zplot)<0.02);
zplot(jj)=0;
cHandle=contourf(x,tArray,zplot);
hold on
if model.k1(1)~=0;plot(model.maxRunupMax_x,model.dates,'w','LineWidth',3);end
colormap(hot)
xlim([x(jx) x(end)])
datetick('y')
c=colorbar;
c.Location='northoutside';
set(gca,'LineWidth',3,'FontSize',24)

disp(['Run-time = ', num2str(toc/60), ' minutes'])
% Program finished
disp(['Program finished at: ',datestr(now)])

return %%%%%%%%% END MAIN PROGRAM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
