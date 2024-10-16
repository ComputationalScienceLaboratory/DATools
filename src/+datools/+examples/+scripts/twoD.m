% run this script to test against some standard time invariant 2D problems
% Problems available
% SineWave, Normal Normal Squared, Normal Gamma (add others)
%
clear;
close all;
clc;

f1 = figure;
f2 = figure;

N = 500; % number of particle/ensemble
problem = 'bimodalGaussian'; % choose the problem
filterName = 'ETPF'; % Run filter of your choice
xf = getEnsembles(problem,N);
observeIndex = 1;
unobservIndex = 2;

% observation operators
HTemp = eye(2);
H = HTemp(observeIndex, :);
y = 2.5;
sigmaObserv = 0.5;
R = sigmaObserv^2 * eye(1);
timeCurrent = 0;

% Observaton model (Gaussian here)
obserrormodel = datools.uncertainty.Gaussian('Covariance', R);
observation = datools.observation.Indexed(2, ...
    'Uncertainty', obserrormodel, ...
    'Indices', observeIndex);
observation.Uncertainty.Mean = y;

modelODE = datools.ODEModel();
modelODE.Y0 = mean(xf,2); % this is the mean of current states
modelODE.TimeSpan = [0, 0]; % this is time invariant

solvermodel = @(f, t, y) ode45(f, t, y);

model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);

filter = datools.filter.ensemble.(filterName)(model, 'InitialEnsemble', xf);

filter.analysis(observation);

xa = filter.Ensemble;

%% plot
% the prior

% the posterior


%% user defined functions
function [xf] = getEnsembles(name,N)
switch name
    case 'SineWave'
        [xf, ~] = datools.toyproblems.twoD.sineWave(N);
    case 'NormalNormalSquared'
        [xf, ~] = datools.toyproblems.twoD.normalNormalsq(N);
    case 'NormalGamma'
        [xf, ~] = datools.toyproblems.twoD.normalGamma(N);
    case 'bimodalGaussian'
        [xf, ~] = datools.toyproblems.twoD.bimodalGaussian(N);
end
end