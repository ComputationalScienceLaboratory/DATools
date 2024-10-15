% this is a standalone script for Lorenz 63
%

clear;
close all;
clc;

% time steps
dt = 0.12;
filterName = 'ETKF';
filterType = 'Ensemble';

% define steps and spinups
spinup = 500;
times = 11 * spinup;

solverModel = @(f, t, y) datools.utils.eDP54(f, t, y);
solverNature = @(f, t, y) datools.utils.eDP54(f, t, y);

% Define ODE (rhs) for truth 
fNature = @(t,x) datools.toyproblems.lorenz63.f(t,x);
natureODE = datools.ODEModel('F', fNature);
natureODE.TimeSpan = [0, dt]; 
nature0 = randn(3, 1);

% Propogate the truth
[tt, yy] = ode45(natureODE.F, [0, 10], nature0);
natureODE.Y0 = yy(end, :).';

% Define ODE (rhs) for the model
fModel = @(t,x) datools.toyproblems.lorenz63.f(t,x);
modelODE = datools.ODEModel('F', fModel);
modelODE.TimeSpan = [0, dt];

% initialize model
model = datools.Model('Solver', solverModel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solverNature, 'ODEModel', natureODE);

% Observation Model
natureToModel = @(x) x;

observeIndicies = 1:1:3;

observSig = 1;

R = (observSig) * speye(numel(observeIndicies));

% Observaton model (Gaussian here)
obsErrorModel = datools.uncertainty.Gaussian('Covariance', R);
observation = datools.observation.Indexed(model.NumVars, ...
    'Uncertainty', obsErrorModel, ...
    'Indices', observeIndicies);

% observation model for the truth
natureObsErrorModel = datools.uncertainty.Gaussian('Covariance', R);
natureObs = datools.observation.Indexed(model.NumVars, ...
    'Uncertainty', natureObsErrorModel, ...
    'Indices', observeIndicies);

modelError = datools.uncertainty.NoUncertainty;

ensembleGenerator = @(N) randn(3, N);

localization = [];
% localization = @(y, H) datools.tapering.bloc.gc(y, r, d, H);
% localization = @(y, H, k) datools.tapering.rloc.gc(y, r, d, H, k);

N = 50;
infl = 1.05;
rej = 0.1;

histVar = 1:1:3;

% define the filter 
filter = datools.filter.ensemble.(filterName)(model, ...
    'InitialEnsemble', ensembleGenerator(N)/10, ...
    'Inflation', infl, ...
    'Parallel', false, ...
    'RankHistogram', histVar, ...
    'Rejuvenation', rej);


filter.MeanEstimate = natureODE.Y0;

mses = zeros(1, times - spinup);
rmse = NaN * ones(1, times-spinup);

for i = 1:times
    % propogate the truth
    nature.evolve();   
    
    % propogate the model
    filter.forecast();

    xt = natureToModel(nature.State);
    y = natureObs.observeWithError(xt);
    observation.Uncertainty.Mean = y;
    
    % Rank histogram (if needed)
    datools.utils.stat.RH(filter, xt, y);

    % analysis
    filter.analysis(observation);
    xa = filter.MeanEstimate;

    if i>spinup
        mses(i - spinup) = mean((xa - xt).^2);
        rmse(i - spinup) = sqrt(mean(mses(1:(i - spinup))));
    end
end





