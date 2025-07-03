% this is a standalone script for Lorenz 96
%
clear;
close all;
clc;

% time steps
dt = 0.05;
F = 8;
numStates = 40; % states
filterName = 'EnKF';
filterType = 'Ensemble';
N = 50; % ensemble
infl = 1.05;
rej = 0.1;

histVar = 1;

% define steps and spinups
spinup = 500;
times = 11 * spinup;

% Time Stepping Methods
solverModel = @(f, t, y) datools.utils.rk4ens(f, t, y, 1);
solverNature = @(f, t, y) datools.utils.rk4ens(f, t, y, 1);

fNature = @(t,x) datools.toyproblems.lorenz96.f(t,x);
natureODE = datools.ODEModel('F', fNature);
natureODE.TimeSpan = [0, dt]; 
natureODE.NumVars = numStates;

nature0 = F*ones(numStates, 1);
nature0(floor(numStates/2)) = 1.001*F; % refer (Lorenz & Emanuel 1998)

% Propogate the truth
[tt, yy] = ode45(natureODE.F, [0, 10], nature0); % can use rk4 method too for speed
natureODE.Y0 = yy(end, :).';

% Define ODE (rhs) for the model
fModel = @(t,x) datools.toyproblems.lorenz96.f(t,x);
modelODE = datools.ODEModel('F', fModel);
modelODE.TimeSpan = [0, dt];
modelODE.NumVars = numStates;

% initialize model
model = datools.Model('Solver', solverModel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solverNature, 'ODEModel', natureODE);

% Observation Model
natureToModel = @(x) x;

observeIndicies = 1:1:numStates;

observSig = 8;

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

ensembleGenerator = @(N) randn(numStates, N); %  this will generate random initialization fo rinitial ensemble



r = 4;
d = @(y, i, j) datools.toyproblems.lorenz96.distanceFunc(0,y,i,j);

localization = [];
% localization = @(y, H) datools.tapering.bloc.gc(y, r, d, H);
% localization = @(y, H, k) datools.tapering.rloc.gc(y, r, d, H, k);

% define the filter 
filter = datools.filter.ensemble.(filterName)(model, ...
    'InitialEnsemble', ensembleGenerator(N)/10, ...
    'Localization', localization, ...
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
        rmse(i - spinup)
    end
    
end
