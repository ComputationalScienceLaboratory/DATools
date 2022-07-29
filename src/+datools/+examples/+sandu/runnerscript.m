% runner script
clear all; close all; clc;
% user defines the filter and everything that are required
% model = Lorenz63
% filtername = 'EnKF';(D)
% filtername = 'ETKF';(D)
% filtername = 'ETPF';
% filtername = 'ETPF2';
% filtername = 'SIR';
% filtername = 'RHF';
% filtername = 'SIS_EnKF';

%
% Lorenz96
% filtername = 'EnKF';
% filtername = 'ETKF';
% filtername = 'LETPF';
% filtername = 'LETKF';
% filtername = 'RHF';

% QG
% filtername = 'EnKF';
% filtername = 'ETKF';
% filtername = 'LETKF';

%% User inputs
options.modelname = 'Lorenz63';

options.filtername = 'EnKF';

options.ensNs = [4, 8, 12, 16, 20, 24, 28];%[25, 50, 75, 100];

options.infs = round(linspace(1, 1.10, 7), 2);    % [1, 1.025, 1.05, 1.075, 1.10];

options.variance = 1; % observation variance

options.observeindicies = 1:1:3;  % observation indices

options.rejs = round(logspace(-1.5, -0.25, 7), 2); % round(2*logspace(-1.5, -0.5, 6), 2)

options.spinups = 500;

options.steps = 11 * options.spinups; % change as deemed fit

options.Dt = 0.12; % 0.12(L63), 0.05(L96), 0.0109(QG)

options.odesolver = 'ode45'; % ode45 , RK4

options.localize = false; % set to true if localization is needed

options.localizationradius = 4; % localization rdius

options.ns = 20;

%plotting parameters
options.rankhistogramplotindex = 1:2:numel(options.ensNs);
options.rmseplotindex = 1:2:numel(options.ensNs);
options.rmseheatmapplotindex = 1:1:numel(options.ensNs);
options.kldivergenceplotindex = 1:1:numel(options.ensNs);

%%
% call the experiment file
datools.examples.sandu.runexperiments2(options);
