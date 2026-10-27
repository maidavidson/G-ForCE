function [Deq,DeqSand,depthAttenuation,overwash,maxRunup_x,maxWaterDepth,Gradient,swashIndex,Ir]=...
    findMeanDissipation_3(x,h,model,wave,tide,seaLevel,startDate,avInterval)
%% Function or program name: findMeanDissipation_3
% Initialise no-overwash
%
%   Inputs:
%       x — input argument (description/units to be verified)
%       h — input argument (description/units to be verified)
%       model — input argument (description/units to be verified)
%       wave — input argument (description/units to be verified)
%       tide — input argument (description/units to be verified)
%       seaLevel — input argument (description/units to be verified)
%       startDate — input argument (description/units to be verified)
%       avInterval — input argument (description/units to be verified)
%
%   Outputs:
%       Deq — output argument (description/units to be verified)
%       DeqSand — output argument (description/units to be verified)
%       depthAttenuation — output argument (description/units to be verified)
%       overwash — output argument (description/units to be verified)
%       maxRunup_x — output argument (description/units to be verified)
%       maxWaterDepth — output argument (description/units to be verified)
%       Gradient — output argument (description/units to be verified)
%       swashIndex — output argument (description/units to be verified)
%       Ir — output argument (description/units to be verified)
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

% Initialise no-overwash
  overwash=0;

% Max and min limits of swash zone
  swashIndexMin=inf;
  swashIndexMax=0;

% Compute start date for the averaging before the model start date
 startDate=startDate-round(avInterval);

% Check that waves exist in order to do computation - estimate going
% forward from the start date if not
 if startDate<wave.dates(1)
  timeVal=wave.dates(1):1/24:wave.dates(1)+round(avInterval); % time-scale for calculations
 else
  timeVal=startDate:1/24:startDate+round(avInterval);
 end

 % Interpolate the wave, tide and sea level conditions over the averaging
 % period
 HYear=interp1(wave.dates,wave.H,timeVal,'linear','extrap'); % corresponding wave height, ...
 TYear=interp1(wave.dates,wave.Tp,timeVal,'linear','extrap'); % period, ...
 thetaYear=interp1(wave.dates,wave.theta,timeVal,'linear','extrap'); % wave angle, ...
 z_tideYear=interp1(tide.dates,tide.z+tide.residual,timeVal,'linear','extrap'); % tide,...
 seaLevelYear=interp1(seaLevel.dates,seaLevel.z,timeVal,'linear','extrap'); % and sea level.

 % Quantities invariant across the averaging loop.
 totalWaterLevelYear=z_tideYear+seaLevelYear;
 hColumn=h(:);
 surfaceD50=model.D50(1);

 % Initialise max run-up and max depth
 maxRunup_x=-inf;
 maxWaterDepth=-inf;
 nPts=length(HYear);
 Gradient=zeros(1,nPts);
 Ir=zeros(1,nPts);

% Main averaging loop here
for ii=1:nPts

% compute the decay of turbulence with depth
  k2=pi./HYear(ii);

% Run wave model in a loop
 % Call wave model
 [waveModel] = WaveModelBJ_swash(HYear(ii),TYear(ii),thetaYear(ii),x,h,totalWaterLevelYear(ii),surfaceD50,model);

% compute swashzone values of dissipation and run-up CSIs
  maxRunupOld=maxRunup_x;
  maxWaterDepthOld=maxWaterDepth;

% record gradient
  Gradient(ii)=waveModel.Gradient;
  Ir(ii)=waveModel.Ir;

% record swashzone limits
  swashIndexMin=min(min(waveModel.swashIndex),swashIndexMin);
  swashIndexMax=max(max(waveModel.swashIndex),swashIndexMax);

% max Run-up and max water depth
  maxRunup_x=max(maxRunupOld,waveModel.runup_x);
  maxWaterDepth=max(maxWaterDepthOld,waveModel.runup_z);

% total water depth
  htot=hColumn+waveModel.eta(:)+totalWaterLevelYear(ii);
  kk=htot<=0;
  htot(kk)=0;
  kk=[];

% depth attenuation of wave generated turbulence
  depthAttenuation=exp(-k2.*max(0,htot));

% Sum swash and attenuated surfzone dissipation
  transportSign=1;%sign(wave.Pmean-P(ii));
  %if Ir>3.3; transportSign=1;end
  Dw= transportSign.*waveModel.Dtot(:).*depthAttenuation;

% Now do thresholding
  if model.layers==2
     DwSand=Dw;
  else
     DwSand=Dw*0;
  end
    [Dcrit]=dissThrehold(TYear(ii), model,1);
    signDw=sign(Dw);
    Dw=(abs(Dw)-Dcrit);
    jj=Dw<0;
    Dw(jj)=0;
    Dw=Dw.*signDw;
% Now sand layer
  if model.layers==2
     [Dcrit]=dissThrehold(TYear(ii), model,2);
     signDwSand=sign(DwSand);
     DwSand=abs(DwSand)-Dcrit;
     jj=DwSand<0;
     DwSand(jj)=0;
     DwSand=DwSand.*signDwSand;
  end

% overwash?
  overwash=max(overwash,waveModel.overtop);

% count finite values of Dw (used for accurate mean calculation)
  countFiniteValues=double(abs(Dw)>0);

% count finite values of sand layer
  countFiniteValues1=double(abs(DwSand)>0);

 % counter for summing finite values of dissipation
 if ii==1
    sumCount=countFiniteValues;
    sumCount1=countFiniteValues1;
 end

% sum values
  if ii==1
     sumD = Dw;
     sumDSand = DwSand;
  else
     sumD = sumD + Dw;
     sumDSand = sumDSand + DwSand;

     sumCount=sumCount+countFiniteValues;
     sumCount1=sumCount1+countFiniteValues1;
  end

end % End averaging loop


% finally calculate the mean
  Deq = sumD  ./sumCount;
  DeqSand = sumDSand  ./sumCount1;

% clear NaN/inf values
  jj=isnan(Deq)==1;
  Deq(jj)=0; clear jj
  jj=isnan(DeqSand)==1;
  DeqSand(jj)=0; clear jj

  jj=isinf(Deq)==1;
  Deq(jj)=0; clear jj
  jj=isinf(DeqSand)==1;
  DeqSand(jj)=0; clear jj

% absolute index of swash zone over the averaging period equitable to Deq
  swashIndex=swashIndexMin:swashIndexMax;

% gradient mean
Gradient=mean(Gradient);
Ir=mean(Ir);

end % End subroutine
