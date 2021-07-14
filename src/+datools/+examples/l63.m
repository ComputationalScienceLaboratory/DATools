clear;
close all;

rng(5);

Deltat =  0.12;

solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 100);
solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 100);

% solvermodel = @(f, t, y) ode45(f, t, y);
% solvernature = @(f, t, y) ode45(f, t, y);

natureODE = otp.lorenz63.presets.Canonical;
natureODE.TimeSpan = [0, Deltat];

nvrs = natureODE.NumVars;

nature0 = randn(nvrs, 1);

[tt, yy] = ode45(natureODE.Rhs.F, [0 10], nature0);
natureODE.Y0 = yy(end, :).';

modelODE = otp.lorenz63.presets.Canonical;
modelODE.TimeSpan = [0, Deltat];


model  = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% naturetomodel = datools.observation.Linear(numel(nature0), 'H', speye(nvrs));
naturetomodel = datools.observation.Indexed(numel(nature0), 'Indices', 1:nvrs);

observeindicies = 1;

nobsvars = numel(observeindicies);

R = (8/1)*speye(nobsvars);

obserrormodel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));
observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indices', observeindicies);

% We make the assumption that there is no model error
modelerror = datools.error.Error;

ensembleGenerator = @(x) randn(nvrs, x);

ensN = 100;
infl = 1.01;
rej = 0.1;

% No localization
localization = [];

meth = datools.statistical.ensemble.SIS_EnKF(model, ...
    'ModelError', modelerror, ...
    'Observation', observation, ...
    'NumEnsemble', ensN, ...
    'EnsembleGenerator', ensembleGenerator, ...
    'Inflation', infl, ...
    'Rejuvenation', rej, ...
    'Localization', localization, ...
    'Parallel', false, ...
    'RIPIterations', 0);

meth.Ensemble = nature.State + meth.Ensemble;

spinup = 0;
times = 1000;

do_enkf = true;

sse = 0;

for i = 1:times
    
    nature.evolve();
    
    if do_enkf
        meth.forecast();
    end
    
    
    % observe
    xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
    y = meth.Observation.observeWithError(model.TimeSpan(1), xt);
    
    % analysis
    %try
        if do_enkf
            meth.analysis(R, y);
        end
    %catch
    %    do_enkf = false;
    %end
    
    xa = meth.BestEstimate;
    
    err = xt - xa;
    
    if i > spinup
        
        sse = sse + sum((xa - xt).^2);
        rmse = sqrt(sse/(i - spinup)/nvrs);
        
    end
    
    if rem(i, 10) == 0
        fprintf('step %d %.5f\n', i, rmse);
    end
    
end

rmse

