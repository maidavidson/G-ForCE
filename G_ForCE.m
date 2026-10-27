function G_ForCE
%% Function or program name: G_ForCE
% Initialise
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

%% Initialise
 clear; close all; clc

%% Setup paths
 addpath([pwd,'/AnalyseForceOutput'])
 addpath([pwd,'/ForCE'])


%% Open new figure maximised
  fig=figure('Name','Forcasting Coastal Evolutuion [ForCE]','NumberTitle','off');
  set(gcf,'units','normalized','outerposition',[0 0 1 1]);
  set(fig,'MenuBar','none')

 %% Open a display window and display welcome text file
 ph = uipanel(fig,'Units','normalized','position', [0.0 0.0 1 0.75],'title',...
    'Help window');
 lbh = uicontrol(ph,'style','listbox','Units','normalized','position',...
    [0.01 0.01 0.99 0.99],'FontSize',14);
 displayText('/help/welcome.txt')

 % Figure position
 pos = [0.05 0.75 0.8 0.25]; % [left bottom width height]
 subplot('Position',pos);
 C = imread('help/startImage.png');
 imshow(C)

%% Create menus
%% Setup
 m = uimenu(fig,'Text','Setup');
 mitem = uimenu(m,'Text','&Open Setup File');
 mitem.Accelerator = 'S';
 mitem.MenuSelectedFcn = @openSetupfun;

 mitem_1 = uimenu(m,'Text','&Prepare spin-up (cold-start) survey.mat file from CCO data');
 mitem_1.Accelerator = 'P';
 mitem_1.MenuSelectedFcn = @prepCCOfun;

 mitem_2 = uimenu(m,'Text','&Prepare initial (hot-start) survey.mat file from CCO data');
 mitem_2.Accelerator = 'I';
 mitem_2.MenuSelectedFcn = @prepCCO4ForCEfun;

 mitem_3 = uimenu(m,'Text','&Generate a synthetic waves file from an existing waves file');
 mitem_3.Accelerator = 'W';
 mitem_3.MenuSelectedFcn = @prepSyntheticWavesfun;

 mitem_4 = uimenu(m,'Text','&Generate a synthetic tides file from an existing tides file');
 mitem_4.Accelerator = 'H';
 mitem_4.MenuSelectedFcn = @prepSyntheticTidesfun;

 mitem_5 = uimenu(m,'Text','&Generate a sea level file');
 mitem_5.Accelerator = 'Z';
 mitem_5.MenuSelectedFcn = @prepSeaLevelfun;

 %% Run
 m1 = uimenu(fig,'Text','&Run');
 mitem1 = uimenu(m1,'Text','&Run ForCE');
 mitem1.Accelerator = 'R';
 mitem1.MenuSelectedFcn = @RunForCEfun;

 %% Calibrate
 % m2 = uimenu(fig,'Text','&Calibrate');
 % mitem2 = uimenu(m2,'Text','&Calibrate ForCE');
 % mitem2.Accelerator = 'C';
 % mitem2.MenuSelectedFcn = @ClalibrateForCEfun;

 %% Analyse
 m3 = uimenu(fig,'Text','&Analyse');
 mitem3_1 = uimenu(m3,'Text','&Plot predicted CSIs');
 mitem3_1.MenuSelectedFcn = @plotCSIfun;

 mitem3_2 = uimenu(m3,'Text','&Plot predicted profile evolution');
 mitem3_2.MenuSelectedFcn = @plotProfileEvolfun;

 mitem3_3 = uimenu(m3,'Text','&Compare with CCO profile');
 mitem3_3.MenuSelectedFcn = @compareWithCCOfun;

 %% Help
 m4 = uimenu(fig,'Text','&Help');
 mitem4_1 = uimenu(m4,'Text','&How to start?');
 mitem4_1.MenuSelectedFcn = @how2Startfun;

 mitem4_2 = uimenu(m4,'Text','&Prepare spin-up (cold-start) survey.mat file from CCO data');
 mitem4_2.MenuSelectedFcn = @HprepCCOfun;

 mitem4_3 = uimenu(m4,'Text','&Prepare initial (hot-start) survey.mat file from CCO data');
 mitem4_3.MenuSelectedFcn = @HprepCCO4ForCEfun;

 mitem4_4 = uimenu(m4,'Text','&Run ForCE');
 mitem4_4.MenuSelectedFcn = @HRunForCEfun;

 mitem4_5 = uimenu(m4,'Text','&Calibrate ForCE');
 mitem4_5.MenuSelectedFcn = @HClalibrateForCEfun;

 mitem4_6 = uimenu(m4,'Text','&Plot predicted CSIs');
 mitem4_6.MenuSelectedFcn = @HplotCSIfun;

 mitem4_7 = uimenu(m4,'Text','&Plot predicted profile evolution');
 mitem4_7.MenuSelectedFcn = @HplotProfileEvolfun;

 mitem4_8 = uimenu(m4,'Text','&Compare with CCO profile');
 mitem4_8.MenuSelectedFcn = @HcompareWithCCOfun;


 %% Menu functions here...
 function openSetupfun(src,event)
 displayText('/help/openSetup.txt')
 %% Get setup file
 [file,path] = uigetfile('*.m');
  if file~=0
   open(fullfile(path,file))
  end
 end

 function RunForCEfun(src,event)
   displayText('/help/run.txt')
   %run('G_ForCE_v1_1.m')
   evalin('base','runG_ForCE');
 end

 % function ClalibrateForCEfun(src,event)
 %  displayText('/help/calibrate.txt')
 %  evalin('base','calibrateForCE');
 % end

 function prepCCOfun(src,event)
  displayText('/help/prepCCO.txt')
  evalin('base','processCCO');
 end

function prepCCO4ForCEfun(src,event)
  displayText('/help/prepCCO4ForCE.txt')
  evalin('base','prepare_CCO_4_ForCE');
end

function plotCSIfun(src,event)
  displayText('/help/plotCSIs.txt')
  evalin('base','plotForCE_CSIs');
end

function plotProfileEvolfun(src,event)
  displayText('/help/profileEvol.txt')
  subplot('Position',pos);
  evalin('base','plotModelProfileEvolution');
end

function compareWithCCOfun(src,event)
  displayText('/help/compareWithCCO.txt')
  evalin('base','compare_ForCE_with_CCO');
end

%% Help functions here....
function how2Startfun(src,event)
  displayText('/help/welcome.txt')
end

 function HRunForCEfun(src,event)
   displayText('/help/run.txt')
 end

 % function HClalibrateForCEfun(src,event)
 %  displayText('/help/calibrate.txt')
 % end

 function HprepCCOfun(src,event)
  displayText('/help/prepCCO.txt')
 end

function HprepCCO4ForCEfun(src,event)
  displayText('/help/prepCCO4ForCE.txt')
end

function HplotCSIfun(src,event)
  displayText('/help/plotCSIs.txt')
end

function HplotProfileEvolfun(src,event)
  displayText('/help/profileEvol.txt')
end

function HcompareWithCCOfun(src,event)
  displayText('/help/compareWithCCO.txt')
end

function prepSyntheticWavesfun(src,event)
  evalin('base','syntheticWaveGenerator');
end

function prepSeaLevelfun(src,event)
  evalin('base','GenerateSL');
end

%% text dispay function
 function displayText(helpFile)
 warning off
 fid = fopen([pwd,helpFile]);
 indic = 1;
 while 1
     tline = fgetl(fid);
     if ~ischar(tline),
         break
     end
     strings{indic}=tline;
     indic = indic + 1;
 end
 fclose(fid);
 strings=textwrap(lbh,strings);
 set(lbh,'string',strings);
 set(lbh,'Value',1);
 set(lbh,'Selected','on');
 set(lbh,'FontSize',18);
 warning on
 end


end
