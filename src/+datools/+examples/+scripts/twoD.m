% run this script to test against some standard time invariant 2D problems
% Problems available
% SineWave, Normal Normal Squared (add others)
%
clear;
close all;
clc;

f1 = figure;
f2 = figure;

N = 200; % number of particle/ensemble
problem = 'SineWave';
filterName = 'EnKF'; % Run filter of your choice
xf = getEnsembles(problem,N);
observeIndex = 1;
unobservIndex = 2;

% observation operators
HTemp = eye(2);
H = HTemp(observeIndex, :);
y = 0.75;
sigmaObserv = 0.25;
R = sigmaObserv^2 * eye(1);
timeCurrent = 0;

% Observaton model (Gaussian here)
obserrormodel = datools.uncertainty.Gaussian('Covariance', R);
observation = datools.observation.Indexed(2, ...
    'Uncertainty', obserrormodel, ...
    'Indices', observeIndex);
observation.Uncertainty.Mean = y;

modelODE = @(x, t) model2D(x, t); % x can be different 2d cases name
solvermodel = @(f, t, y) ode45(f, t, y);

model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);

filter = datools.filter.ensemble.(filterName)([], 'InitialEnsemble', xf);

filter.analysis(observation);

xa = filter.MeanEstimate;

%% plot
% the prior

% the posterior


%% user defined functions
function [xf] = getEnsembles(name,N)
switch name
    case 'SineWave'
        xf = datools.examples.twoDProblems.sineWave(N);
    case 'Normal Normal Squared'
        xf = datools.examples.twoDProblems.normalNormalsq(N);
end
end