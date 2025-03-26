clear;
close all;
clc;

rng(1729);

% time steps
dt = 0.12;
filterName = 'ThreeDVar';
filterType = 'Particle';

Deltat = 0.12;

% solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 100);
% solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 100);
%
% solvermodel = @(f, t, y) ode45(f, t, y);
% solvernature = @(f, t, y) ode45(f, t, y);

solverModel = @(f, t, y) datools.utils.eDP54(f, t, y);
solverNature = @(f, t, y) datools.utils.eDP54(f, t, y);

% Define ODE for truth
otpNature = otp.lorenz63.presets.Canonical;
% natureODE = otp.lorenz63.presets.Canonical;
natureODE = datools.ODEModel('OTPObject', otpNature);
nature0 = randn(natureODE.NumVars, 1);
natureODE.TimeSpan = [0, dt];

% Define ODE for the model
otpModel = otp.lorenz63.presets.Canonical;
% modelODE = otp.lorenz63.presets.Canonical;
modelODE = datools.ODEModel('OTPObject', otpModel);
modelODE.TimeSpan = [0, dt];

% natureODE = otp.lorenz63.presets.Canonical;
% natureODE.TimeSpan = [0, Deltat];
%
nvrs = natureODE.NumVars;
%
% nature0 = randn(nvrs, 1);
%
% [tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
% natureODE.Y0 = yy(end, :).';
%
% modelODE = otp.lorenz63.presets.Canonical;
% modelODE.TimeSpan = [0, Deltat];

% Propogate the truth
% [tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
[tt, yy] = ode45(natureODE.F, [0, 10], nature0);
natureODE.Y0 = yy(end, :).';

% We make the assumption that there is no model error
modelUncertainty = datools.uncertainty.NoUncertainty;

% initialize model
model = datools.Model('Solver', solverModel, 'ODEModel', modelODE, 'Uncertainty', modelUncertainty);
nature = datools.Model('Solver', solverNature, 'ODEModel', natureODE, 'Uncertainty', modelUncertainty);

% model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
% nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% naturetomodel = datools.observation.Linear(numel(nature0), 'H', speye(nvrs));
% naturetomodel = datools.observation.Indexed(numel(nature0), 'Indices', 1:nvrs);
natureToModel = @(x) x;

% observeindicies = 1:3;
% observe these variables
observeIndicies = 1:1:natureODE.NumVars;

nobsvars = numel(observeIndicies);

R = (8 / 1) * speye(nobsvars);
dR = decomposition(R, 'chol');

% obserrormodel = datools.uncertainty.Gaussian('Covariance', R);
% observation = datools.observation.Indexed(model.NumVars, ...
%     'Uncertainty', obserrormodel, ...
%     'Indices', observeindicies);

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

% We make the assumption that there is no model error

% No localization
localization = [];

B = [0.86, 0.86, -0.02; ...
    0.86, 1.1149, -0.01; ...
    -0.02, -0.01, 1.02];

B = 8 * B;

modelerror = datools.uncertainty.Gaussian('Covariance', B/2);


meth = datools.filter.variational.(filterName)(model, ...
    'CovarianceEstimate', [], ...
    'InitialState', nature.State, ...
    'BackgroundCovariance', B, ...
    'OptimizationType', 'lbfgs');

% meth = datools.filter.gaussian.UKF(model, ...
%     'ModelError', modelerror, ...
%     'Observation', observation, ...
%     'InitialState', nature.State, ...
%     'InitialCovariance', B, ...
%     'Alpha', 1e-1);

spinup = 500;
times = 11 * spinup;

do_filter = true;

sse = 0;
rmse = 0;

for i = 1:times

    nature.evolve();

    if do_filter
        meth.forecast();
    end


    % observe
    % xt = naturetomodel.observeWithoutError(nature.State);
    % y = meth.Observation.observeWithError(xt);
    xt = natureToModel(nature.State);
    y = natureObs.observeWithError(xt);
    observation.Uncertainty.Mean = y;
    

    % analysis
    %try
    if do_filter
        meth.analysis(observation);
    end
    %catch
    %    do_enkf = false;
    %end

    xa = meth.MeanEstimate;

    err = xt - xa;

    if i > spinup

        sse = sse + sum((xa - xt).^2);
        rmse = sqrt(sse/(i - spinup)/nvrs);

    end

    if rem(i, 10) == 0
        fprintf('step = %d rmse = %.5f\n', i, rmse);
    end

end
