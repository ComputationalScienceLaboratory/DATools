% runner script
clear all; close all; clc;
% user defines the filter and everything that are required
% model = Lorenz63
% filtername = 'EnKF';
% filtername = 'ETKF';
% filtername = 'ETPF';
% filtername = 'ETPF2';
% filtername = 'SIR';
% filtername = 'RHF';
%
% Lorenz96
% filtername = 'EnKF';
% filtername = 'ETKF';
% filtername = 'LETPF';
% filtername = 'SIR';(X)
% filtername = 'LETKF';
% filtername = 'RHF';

% QG
% filtername = 'EnKF';
% filtername = 'ETKF';
% filtername = 'LETKF';

%% User inputs
options.modelname = 'Lorenz63';

options.filtername = 'RHF';

options.ensNs = [10, 20, 30, 40, 50, 100];%[25, 50, 75, 100];

options.infs = [1, 1.025, 1.05, 1.075, 1.10];

options.variance = 1; % observation variance

options.rejs = round(2*logspace(-1.5, -0.5, 6), 2);

options.spinups = 500;

options.steps = 11 * options.spinups; % change as deemed fit

options.Dt = 0.12; % 0.12(L63), 0.05(L96), 0.0109(QG)

options.odesolver = 'ode45'; % ode45 , RK4

options.localize = false;

options.ns = 20;

%%
% call the experiment file
datools.examples.sandu.runexperiments2(options);
