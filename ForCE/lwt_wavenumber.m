%%%%%
function k=lwt_wavenumber(f,h,del_k, k_deep);
%% Function or program name: lwt_wavenumber
% LWT_WAVENUMBER G_ForCE support routine; see function signature for inputs and outputs.
%
%   Inputs:
%       f — input argument (description/units to be verified)
%       h — input argument (description/units to be verified)
%       del — input argument (description/units to be verified)
%       kguess — input argument (description/units to be verified)
%
%   Outputs:
%       k — output argument (description/units to be verified)
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

%LWT_WAVENUMBER G_ForCE support routine; see function signature for inputs and outputs.

% This code was written by Mark Davidson and amended for application to gravel and mixed beaches as part of the NERC Gravel Beach Project.
% function k=lwt_wavenumber(f,h,del, kguess);
%
% linear gravity wave dispersion relationship is iterated to calculate wave number at particular frequencies
%
% REMEMBER if you are plotting a dispersion curve, simply calculate f=sqrt(g k tanh(k h))/(2*pi)
%
%  k=lwt_wavenumber(f,h,del, kguess);
%
% input
%    f is frequency
%    h is water depth
%    del is tolerance (in terms of k) to iterate to, a good choice is ko/100
%    kguess is an initial guess
%
% output
%    k is wavenumber = 2 Pi/L

nf = max(size(f));
g = 9.81;

for i = 1:nf

% start out with deep water approx
    ko = (4*pi*pi*f(i)^2)/g;
    if (kguess<0)
        kguess = ko;
    end

% check for asymptotic conditions
% deep water
    if (h*ko/(2*pi) > 0.5)
        k(i) = ko;
%        fprintf('deep water\n')

% shallow water
    elseif (h*ko/(2*pi) < 0.005)
        k(i) = 2*pi*f/sqrt(g*h);
%        fprintf('shallow water\n')

% otherwise, start iteration
    else
%        fprintf('intermediate water\n')
        del_k = 1e35;

        while (abs(del_k) > del)
            k(i) = ko/tanh(kguess*h);
            del_k = kguess-k(i);
            kguess = 0.5*(k(i)+kguess);
        end
    end
end
end

%%%%%%%%%
