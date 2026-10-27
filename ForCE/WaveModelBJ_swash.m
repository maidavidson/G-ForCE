function waveModel = WaveModelBJ_swash(Hs, T, theta0, x, h, z_tide, D50, model)
%% Function or program name: WaveModelBJ_swash
% WAVEMODELBJ_SWASH Cross-shore wave transformation including swash dynamics.
%
%   Inputs:
%       Hs — input argument (description/units to be verified)
%       T — input argument (description/units to be verified)
%       theta0 — input argument (description/units to be verified)
%       x — input argument (description/units to be verified)
%       h — input argument (description/units to be verified)
%       z_tide — input argument (description/units to be verified)
%       D50 — input argument (description/units to be verified)
%       model — input argument (description/units to be verified)
%
%   Outputs:
%       waveModel — output argument (description/units to be verified)
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
%WAVEMODELBJ_SWASH Cross-shore wave transformation including swash dynamics.
%
% Computes wave height, setup, wave/roller energy, dissipation, reflection,
% run-up and related parameters using a modified Battjes-Janssen model.
% Includes roller evolution, infiltration, bed friction and swash-zone
% processes.
%
% This version focuses on readability and low-risk efficiency improvements.
% Validate against the original before using results scientifically.

%% Constants and model metadata
g = 9.81;
rho = 1027;

model.waveModel = 'modified Battjes-Jansen wave model to include swashzone dynamics';

% Optional model coefficients and defaults.
alpha = getModelValue(model, 'alpha', 0.75);
beta  = getModelValue(model, 'beta', 0.1);
gamma = getModelValue(model, 'gamma', []);
hmin  = getModelValue(model, 'hmin', 0.01);
ksw   = getModelValue(model, 'ksw', 0.025);

%% Input preparation
Hrms0 = Hs ./ sqrt(2);
Nx = numel(x);
dx = diff(x);

% Limit extreme incident angles to avoid downstream numerical problems.
if abs(theta0) > 85
    theta0 = 85 * sign(theta0);
end
theta0 = deg2rad(theta0);

%% Deep-water and offshore boundary conditions
L_deep = g * T^2 / (2*pi);
k_deep = 2*pi / L_deep;
del_L = 0.001;
del_k = 1/L_deep - 1/(L_deep + del_L);
tol = 1e-3;

c_deep  = L_deep / T;
cg_deep = 0.5 * c_deep;
h_tot0 = h(1) + z_tide;

k0 = lwt_wavenumber(1/T, h_tot0, del_k, k_deep);
L0 = 2*pi / k0;
c0 = L0 / T;

sinThetaDeep = c_deep * sin(theta0) / c0;
if abs(sinThetaDeep) < 1
    theta_deep = asin(sinThetaDeep);
else
    theta_deep = theta0;
end

cg0 = (0.5 + k0*h_tot0/sinh(2*k0*h_tot0)) * c0;
Hrms_deep = sqrt(Hrms0^2 * cg0*cos(theta0) / ...
                 (cg_deep*cos(theta_deep)));
s_deep = Hrms_deep / L_deep;

if isempty(gamma)
    gamma = 0.5 + 0.4*tanh(33*s_deep);
end

%% Preallocate arrays
wi = infiltrationVelocity(D50, [], [], []);

c           = zeros(Nx, 1);
cg          = zeros(Nx, 1);
k           = zeros(Nx, 1);
Hrms        = zeros(Nx, 1);
TWL         = z_tide * ones(Nx, 1);
Hroller     = zeros(Nx, 1);
Htot        = zeros(Nx, 1);
E           = zeros(Nx, 1);
Er          = zeros(Nx, 1);
P           = zeros(Nx, 1);
Dw          = zeros(Nx, 1);
Dtot        = zeros(Nx, 1);
Dbf         = zeros(Nx, 1);
DF1         = zeros(Nx, 1);
DF2         = zeros(Nx, 1);
Dr          = zeros(Nx, 1);
Di          = zeros(Nx, 1);
eta         = zeros(Nx, 1);
swashWeight = zeros(Nx, 1);
h_tot       = zeros(Nx, 1);
h_tot2      = zeros(Nx, 1);
Qb          = zeros(Nx, 1);
theta       = zeros(Nx, 1);
n           = zeros(Nx, 1);
Sxx         = zeros(Nx, 1);

QbLast = Qb(1);
hmax = inf;
R2 = 0;

%% Initialise seaward boundary
Hrms(1) = Hrms0;
Htot(1) = Hrms0;
E(1) = 0.125 * rho * g * Hrms(1)^2;
Er(1) = 0;
eta(1) = 0;

h_tot(1) = h(1) + z_tide + eta(1);
TWL(1) = z_tide;
h_tot2(1) = h_tot(1);

k(1) = k0;
c(1) = c0;
cg(1) = cg0;
theta(1) = theta0;

n(1) = 0.5 * (1 + 2*k(1)*h_tot(1)/sinh(2*k(1)*h_tot(1)));
Sxx(1) = ((cg(1)/c(1)) * (1 + cos(theta(1))^2) - 0.5) * E(1) ...
         + 2*Er(1)*cos(theta(1))^2;

waveModel.overtop = 0;
R = 0;
ilast = 1;
swashFlag = 0;
iCrest = find(-h == max(-h), 1, 'first');

%% Identify swash zone and beach gradient
tanBeta = gradient(-h, x);

minWaterDepth = -0.35 * Hs(1) + z_tide;
maxWaterDepth =  0.5 * Hs(1) + z_tide;

j1 = find(-h >= minWaterDepth, 1, 'first');
j2 = find(-h >= maxWaterDepth, 1, 'first');

if isempty(j1)
    x0 = x(end);
    Gradient = 0;
    j1 = numel(x) - 1;
    j2 = j1 + 1;
else
    if isempty(j2)
        j2 = j1 + 1;
    end
    jj = j1:j2;
    x0 = x(jj(1));
    Gradient = mean(tanBeta(jj));
    Gradient = max(min(Gradient, 0.7), 0.01);
end

swashWidth = (maxWaterDepth - minWaterDepth) / Gradient;
kswash = 2*pi / swashWidth;
swashWeight = tanh(kswash .* max(x - x0, 0));

%% Initial reflection coefficient and Iribarren number
Ir = IrribarenNo(Hs(1), T, Gradient);
R0 = reflectionCoefficientSA(Ir, 1, 5);

%% Cross-shore propagation loop
for i = 1:(Nx - 1)

    % Flag overtopping once the profile crest has been passed.
    if i > iCrest
        waveModel.overtop = 1;
    end

    % Water depth at the next point, assuming slowly varying setup.
    if -h(i+1) > z_tide || (i > iCrest && z_tide < max(-h))
        h_tot(i+1) = eta(i);
        swashFlag = 1;
    else
        h_tot(i+1) = h(i+1) + z_tide + eta(i);
    end


    % Wavenumber and effective water depth.
    k(i+1) = lwt_wavenumber(1/T, max(hmin, h_tot(i)), del_k, k(i));
    h_tot(i+1) = max(0, h_tot(i+1));

    wt = tanh(k(i+1) * h_tot(i+1));
    h_tot2(i+1) = h_tot(i+1) + (1 - wt) * (Hrms(i)/2);
    TWL(i+1) = -h(i+1) + h_tot2(i+1);

    n(i+1) = 0.5 * (1 + 2*k(i+1)*h_tot(i+1) / ...
                    sinh(2*k(i+1)*h_tot(i+1)));
    n(i+1) = max(0, min(1, n(i+1)));

    c(i+1) = (g*T)/(2*pi) * tanh(k(i+1)*h_tot2(i+1));
    cg(i+1) = max(0, n(i+1)*c(i+1));
    theta(i+1) = asin(c(i+1)*sin(theta(i))/c(i));

    %% Breaking-wave fraction
    Hm = (0.88/k(i+1)) * tanh(gamma*k(i+1)*h_tot(i)/0.88);

    if i == 1
        Qb1 = 0.005;
    else
        Qb1 = Qb(i-1);
    end

    waveratio = -(Hrms(i)/Hm)^2;

    if 1 - Qb1 < tol
        Qb(i) = 0.999;
    else
        nIter = 0;
        while true
            nIter = nIter + 1;
            Qb2 = Qb1;

            fQb = Qb2 - exp((1 - Qb2)/waveratio);
            dfQbdQb = 1 + exp((1 - Qb2)/waveratio)/waveratio;
            Qb1 = Qb2 - fQb/dfQbdQb;

            if Qb1 < 0
                Qb1 = 0;
            end
            if Qb1 > 1
                Qb1 = 0.999;
            end
            if nIter > 10
                Qb1 = QbLast;
            end

            if abs(Qb2 - Qb1) < tol
                Qb(i) = Qb1;
                break
            end
        end
        QbLast = Qb(i);
    end

    %% Breaking, bore-collapse and infiltration dissipation
    Dw1 = (1 - swashWeight(i)) * 0.25 * alpha * Qb(i) * rho * ...
          (1/T) * g * Hm^2;

    Dw2 = swashWeight(i) * 2*sin(beta) * g * E(i) / c(i);
    if ~isfinite(Dw2)
        Dw2 = 0;
    end

    Di(i) = rho * g * mean(h_tot2(i:i+1), 'omitnan') * wi * swashWeight(i);

    %% Bottom-friction dissipation
    sinhKh = sinh(k(i) * max(hmin, h_tot2(i)));
    Adelta = Hrms(i) / (2*sinhKh);
    fw = exp(-6 + 5.2*(Adelta/ksw)^(-0.19));
    fw = min(fw, 0.3);

    DF1(i) = (1 - swashWeight(i)) * rho * fw * ...
             ((2*pi*Hrms(i))/(T*sinhKh))^3 / (12*pi);

    DF2(i) = swashWeight(i) * 0.001 * rho * c(i)^3 * ...
             (1 / cos(atan(tanBeta(i))));

    Dbf(i) = DF1(i) + DF2(i);

    Dsurf_i = Dw1 + DF1(i);
    Dswash_i = Dw2 + Di(i) + DF2(i);
    Dtot(i) = Dsurf_i + Dswash_i;

    %% Roller dissipation
    Dr(i) = 2*sin(beta) * g * Er(i) / c(i);
    if ~isfinite(Dr(i))
        Dr(i) = 0;
    end

    %% Reflection and reflected-energy cutoff
    % Accumulate the trapezoidal integral incrementally instead of
    % re-integrating the entire offshore-to-i interval at every point.
    % Note: this changes floating-point summation order; validate against
    % the original before accepting the optimisation.
    if i > 2
        if i == 3
            EtotRunning = trapz(x(1:2)', DF2(1:2)' + Di(1:2)');
        end
        EtotRunning = EtotRunning + (x(i)-x(i-1)) * ...
            ((DF2(i-1)+Di(i-1)) + (DF2(i)+Di(i))) / 2;
        Etot = EtotRunning;
        Rsquared = R0^2 * (1 - Etot/E(1));
        Rsquared = max(0, Rsquared);
        R = max(0.01, sqrt(Rsquared));
        Hstop = R * Hrms(1);
    else
        Hstop = 0;
    end

    %% Wave energy balance
    E(i+1) = (E(i)*cg(i)*cos(theta(i)) - Dtot(i)*dx(i)) / ...
             (cg(i+1)*cos(theta(i+1)));
    E(i+1) = max(0, E(i+1));

    %% Roller energy balance
    if E(i) > 0
        dEr = Dw1 - Dr(i);
    else
        dEr = 0;
    end

    Er(i+1) = (2*Er(i)*c(i)*cos(theta(i)) + dEr*dx(i)) / ...
              (2*c(i+1)*cos(theta(i+1)));
    Er(i+1) = max(Er(i+1), 0);
    Hroller(i+1) = sqrt(8*Er(i+1)/(rho*g));

    %% Wave power, wave height and setup
    if i == 1
        P(i) = E(i) * cg(i);
    end
    P(i+1) = (E(i+1) + Er(i+1)) * cg(i+1);

    Hrms(i+1) = sqrt(8*E(i+1)/(rho*g));
    Htot(i+1) = hypot(Hrms(i+1), Hroller(i+1));

    Sxx(i+1) = ((cg(i+1)/c(i+1)) * ...
                (1 + cos(theta(i+1))^2) - 0.5) * E(i+1) ...
                + 2*Er(i+1)*cos(theta(i+1))^2;

    dSxxdx = (Sxx(i+1) - Sxx(i)) / dx(i);

    if swashFlag == 1
        deta = 0;
    else
        deta = -dSxxdx * dx(i) / ...
               (rho*g*max(0.1, mean(h_tot(i:i+1), 'omitnan')));
    end

    if isnan(deta)
        deta = 0;
    end
    eta(i+1) = eta(i) + deta;

    %% Stop at reflected-energy cutoff or run-up limit
    if Hrms(i+1) <= Hstop
        Dtot(i) = 0;
        TWL(i+1:end) = -h(i+1:end);
        break
    end

    ilast = i + 1;
end

%% Complete initial values and convert units
Dw(1) = Dw(2);
Dbf(1) = Dbf(2);

theta = rad2deg(theta);
Hs = Hrms * sqrt(2);
Hroller = Hroller * sqrt(2);
Htot = Htot * sqrt(2);

% Shift dissipation back onto the model grid.
Dtot(1) = Dtot(2);
dxFirst = x(2) - x(1);
Dtot = interp1(x(:) + dxFirst/2, Dtot(:), x, 'linear', 'extrap');

%% Assemble output structure
waveModel.maxWL = TWL';
waveModel.alpha = alpha;
waveModel.beta = beta;
waveModel.gamma = gamma;
waveModel.T = T;
waveModel.D50 = D50;
waveModel.x = x';
waveModel.h = h';
waveModel.z_tide = z_tide;
waveModel.Hs = Hs;
waveModel.Hroller = Hroller;
waveModel.Htot = Htot;
waveModel.theta = theta;
waveModel.eta = eta;
waveModel.Power = P;
waveModel.cg = cg;
waveModel.c = c;
waveModel.k = k;
waveModel.Sxx = Sxx';
waveModel.n = n';
waveModel.Dw = Dw;
waveModel.Dtot = Dtot';
waveModel.Dbf = Dbf;
waveModel.Di = Di;
waveModel.Droller = Dr;
waveModel.Qb = Qb;

% double check R2 value obtained
Gradient= sum(tanBeta(:).*Dw(:))./sum(Dw);
Gradient = min(max(Gradient, 0.01), 0.4);
waveModel.Gradient = Gradient;
if model.swash==1
    [R2]=Poate(Hs(1),T,Gradient, model.swashFact);
else
    [R2,~,~]=StockdonHolman(Hs(1),T,Gradient,model.swashFact);
end
hmax=z_tide+R2;

waveModel.R0 = R0;
waveModel.R = R;
waveModel.Hstop = Hstop;
waveModel.H0 = Hs(1);

iNeeded = min(Nx, ilast + 1);
waveModel.closureError = Hrms(iNeeded) - waveModel.R*Hrms(1);
waveModel.Ir = Ir;

if waveModel.overtop == 1
    waveModel.runup_z = max(TWL);
    waveModel.runup_x = x(iNeeded);
else
    waveModel.runup_z = -h(ilast + 1);
    waveModel.runup_x = x(ilast + 1);
end

if isempty(waveModel.runup_x)
    waveModel.runup_x = x(iNeeded);
end

jj = find(-h > z_tide, 1, 'first');
i0 = jj - 1;
i1 = jj;
xo = [-h(i0), -h(i1)];
yo = [x(i0), x(i1)];
waveModel.shoreline_x = myInterp(xo, yo, z_tide);

waveModel.R2 = waveModel.runup_z - z_tide;
waveModel.R2p = R2;

if isreal(eta)
    waveModel.setup = interp1(x, eta, waveModel.shoreline_x);
else
    waveModel.setup = 0;
end

swashMin = z_tide - waveModel.R2 + waveModel.setup;
if isempty(swashMin)
    i0 = find(-h == max(-h), 1, 'first');
else
    i0 = find(-h >= swashMin, 1, 'first');
end


if waveModel.overtop==1
    i1=min(ilast+1,length(x));
else
    i1=min(ilast+1,length(x));
    if isempty(i1);i1=1;end
    i2=find(-h>=hmax,1,'first');
    if isempty(i2);i2=1;end
    i1=max([i1,i2]);
end

waveModel.swashIndex = (i0:i1)';
waveModel.swashWidth = swashWidth;
waveModel.kswash = kswash;
waveModel.wi = wi;

end


%% Local helper functions

function value = getModelValue(model, fieldName, defaultValue)
%GETMODELVALUE Return a model field or its default.
if isfield(model, fieldName)
    value = model.(fieldName);
else
    value = defaultValue;
end
end


function k = lwt_wavenumber(f, h, del, kguess)
%LWT_WAVENUMBER Solve the linear-wave dispersion relation for wavenumber.
g = 9.81;
nf = max(size(f));

for i = 1:nf
    ko = (4*pi^2*f(i)^2) / g;

    if kguess < 0
        kguess = ko;
    end

    if h*ko/(2*pi) > 0.5
        k(i) = ko;
    elseif h*ko/(2*pi) < 0.005
        k(i) = 2*pi*f / sqrt(g*h);
    else
        del_k = inf;
        while abs(del_k) > del
            k(i) = ko / tanh(kguess*h);
            del_k = kguess - k(i);
            kguess = 0.5 * (k(i) + kguess);
        end
    end
end
end


function Ir = IrribarenNo(Ho, Tp, tanb)
%IRRIBARENNO Calculate the Iribarren number.
if Ho == 0
    Ir = inf;
    return
end
Lo = LoLinear(Tp/1.1);
Ir = tanb / sqrt(Ho/Lo);
end


function R = reflectionCoefficientSA(Ir, a, b)
%REFLECTIONCOEFFICIENTSA Calculate the reflection coefficient.
if isinf(Ir)
    R = a;
    return
end
Ir2 = Ir.^2;
R = (a*Ir2) / (b + Ir2);
R = min(max(R, 0), 1);
end


function wi = infiltrationVelocity(D50, C, hg, poreSpace)
%INFILTRATIONVELOCITY Estimate infiltration velocity [m/s].
if isempty(C)
    C = 0.025;
end
if isempty(hg)
    if D50 < 0.002
        hg = 0.5;
    else
        hg = 1.5;
    end
end
if isempty(poreSpace)
    if D50 < 0.002
        poreSpace = 0.4;
    else
        poreSpace = 0.35;
    end
end
K = C * D50.^2;
q = K * hg;
wi = q / poreSpace;
end


function Lo = LoLinear(T)
%LOLINEAR Calculate deep-water wavelength [m].
g = 9.81;
Lo = g*T.^2 / (2*pi);
end


function R = Poate(Hs, Tp, tanb, fact)
%POATE Run-up parameterisation for gravel beaches.
R = 0.33 .* fact .* sqrt(abs(tanb)) .* Tp .* Hs;
end


% NOTE: Both setup and swash are scaled by swashFact, matching the
% standalone StockdonHolman.m implementation.
