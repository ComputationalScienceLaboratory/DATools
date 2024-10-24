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
options.modelname = 'QG';

options.filtername = 'ETKF';

%options.ensNs = [8, 12, 16, 20, 24, 28, 32];%[25, 50, 75, 100];
options.ensNs = [16, 32, 48, 64, 80, 96, 112];
%options.ensNs = [32, 48, 64, 80, 96, 112, 128];
%options.ensNs = [32, 64, 128, 256, 512, 1024, 2048];


options.infs = round(linspace(1, 1.10, 7), 2);    % [1, 1.025, 1.05, 1.075, 1.10];

options.variance = 1; % observation variance

options.observeindicies = linspace(0,8001,150);  % observation indices

options.rejs = round(logspace(-1.5, -0.25, 7), 2); % round(2*logspace(-1.5, -0.5, 6), 2)

options.spinups = 50;

options.steps = 11 * options.spinups; % change as deemed fit

options.Dt = 0.0109; % 0.12(L63), 0.05(L96), 0.0109(QG)

options.odesolver = 'ode45'; % ode45 , RK4

options.localize = false; % set to true if localization is needed

options.localizationradius = 4; % localization radius

options.ns = 10;

options.histvar = 1:1:1;
options.RHmeasure = 'Truth';

%plotting parameters
options.rankhistogramplotindex = 1:2:numel(options.ensNs);
options.rmseplotindex = 1:2:numel(options.ensNs);
options.rmseheatmapplotindex = 1:2:numel(options.ensNs);
options.kldivergenceplotindex = 1:2:numel(options.ensNs);

%%
% call the experiment file
datools.examples.sandu.runexperiments2(options);
