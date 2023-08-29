clear;
close all;

rng(5);

Deltat = 0.12;

solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 100);
solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 100);

solvermodel = @(f, t, y) ode45(f, t, y);
solvernature = @(f, t, y) ode45(f, t, y);

natureODE = otp.lorenz63.presets.Canonical;
natureODE.TimeSpan = [0, Deltat];

nvrs = natureODE.NumVars;

nature0 = randn(nvrs, 1);

[tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
natureODE.Y0 = yy(end, :).';

modelODE = otp.lorenz63.presets.Canonical;
modelODE.TimeSpan = [0, Deltat];


model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% naturetomodel = datools.observation.Linear(numel(nature0), 'H', speye(nvrs));
naturetomodel = datools.observation.Indexed(numel(nature0), 'Indices', 1:nvrs);

observeindicies = 1:3;

nobsvars = numel(observeindicies);

R = (8 / 1) * speye(nobsvars);
dR = decomposition(R, 'chol');

obserrormodel = datools.uncertainty.Gaussian('Covariance', R);
observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indices', observeindicies);

% We make the assumption that there is no model error

% No localization
localization = [];

B = [0.86, 0.86, -0.02; ...
    0.86, 1.1149, -0.01; ...
    -0.02, -0.01, 1.02];

B = 8*B;

modelerror = datools.uncertainty.Gaussian('Covariance', B/2);


% meth = datools.variational.ThreeDVar(model, ...
%     'ModelError', modelerror, ...
%     'Observation', observation, ...
%     'InitialState', nature.State, ...
%     'BackgroundCovariance', B, ...
%     'OptimizationType', 'lbfgs');

meth = datools.gaussian.UKF(model, ...
    'ModelError', modelerror, ...
    'Observation', observation, ...
    'InitialState', nature.State, ...
    'InitialCovariance', B, ...
    'Alpha', 1e-1);

spinup = 500;
times = 5500;

do_filter = true;

sse = 0;
rmse = 0;

for i = 1:times

    nature.evolve();

    if do_filter
        meth.forecast();
    end


    % observe
    xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
    y = meth.Observation.observeWithError(model.TimeSpan(1), xt);

    % analysis
    %try
    if do_filter
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
        fprintf('step = %d rmse = %.5f\n', i, rmse);
    end

end

rmse
