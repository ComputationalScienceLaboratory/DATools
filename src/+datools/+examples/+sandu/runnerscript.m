% runner script
clear all; close all; clc;
% user defines the filter and everything that are required
% model = Lorenz63
% filtername = 'EnKF';
% filtername = 'ETKF';
% filtername = 'ETPF';
% filtername = 'SIR';
% filtername = 'RHF';
%
% Lorenz96
% filtername = 'EnKF';
% filtername = 'ETKF';
% filtername = 'ETPF';
% filtername = 'SIR';
% filtername = 'LETKF';
% filtername = 'RHF';

%% User inputs
user.modelname = 'Lorenz63';

user.filtername = 'SIS_EnKF';

%user.ensNs = [25, 50, 75, 100];
user.ensNs = [1000];

%user.infs = [1.01, 1.02, 1.05];
user.infs = [1.01];

user.variance = 1; % observation variance

%user.rejs = round(2*logspace(-2, -1, 4), 2);
user.rejs = [0.1];

user.spinups = 2;

user.steps = 11 * user.spinups; % change as deemed fit

user.Dt = 0.12; % 0.12(L63), 0.05(L96)

user.odesolver = 'ode45'; % ode45 , RK4

user.localize = false;

%%
% call the experiment file
datools.examples.sandu.runexperiments(user);
