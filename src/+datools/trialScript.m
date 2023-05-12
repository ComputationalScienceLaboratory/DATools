% trial runner script
clc;
clear all;
close all;

%% Preliminaries

% integration time-step
Delta_t = 0.05;

% filtername
filtername = 'LETKF';

% model
modelName = 'lorenz96'; 

% define steps and spinups
spinup = 500;
steps = 11 * spinup;

% Set number of samples for different initialization
numSample = 1;

% use localization(or not) and radius
localize = true;
r = 4;

% observation error covariance
variance = 8;

modelError = datools.error.Gaussian; % we ignore model error (for now)
addErrorToModel = false; % default is false (for now)


%% define the ODE Model
odeOTP = otp.(modelName).presets.Canonical;
odeModel = datools.ODEModel('OTPObject', odeOTP);
odeModel.TimeSpan = [0, Delta_t];

% define ODE solvers (time integrators) for truth and model
solverModel = @(f, t, y) ode45(f, t, y);
solverTruth = @(f, t, y) ode45(f, t, y);

%define the truth and model ODE
model = datools.Model('Solver', solverModel, 'ODEModel', odeModel, 'SynthError', modelError,...
    'AddError', addErrorToModel);
truth = datools.Model('Solver', solverTruth, 'ODEModel', odeModel, 'SynthError', modelError,...
    'AddError', addErrorToModel);


%% define observation and model error

% observe these variables
observeIndicies = 1:10:40; % change accordingly

numStateObserved = numel(observeIndicies);

R = variance * speye(numStateObserved);

% Observaton error model (Gaussian here)
obsErrorModel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));

% Observation Operator
observationOperator = datools.observation.operator.IndexedObservation(model.NumVars, ...
    'ErrorModel', obsErrorModel, ...
    'Indices', observeIndicies);

% Observation
observation = datools.observation.StateObservation('Covariance', R, 'Noise', obsErrorModel, ...
    'NumObs', numStateObserved);

% Observation model to map from state space to observation space
stateToObservation = datools.observation.operator.LinearObservation(truth.NumVars, ...
    'H', speye(truth.NumVars));

% create random ensemble generator
ensembleGenerator = @(N) randn(model.NumVars, N);


%% define ensemble and inflation
ensNs = [5, 15, 25, 50];
infs = [1.01, 1.02, 1.05, 1.10];
rejs = round(logspace(-1.5, -0.25, 7), 2);

rmses = inf * ones(numel(ensNs), numel(infs));

runsLeft = find(rmses == inf);

for runn = runsLeft.'
    [ensNi, infi] = ind2sub([numel(ensNs), numel(infs)], runn);

    %numSample = 1;
    sE = zeros(numSample, 1);

    inflationAll = infs(infi);

    ensN = ensNs(ensNi);

    for sample = 1:numSample
        % Set rng for standard experiments
        rng(17+sample-1);

        fprintf('N: %d | inf: %.3f | sample = %d\n', ensNs(ensNi), infs(infi), sample);

        inflation = inflationAll;

        if (localize)
            d = @(t, y, i, j) modelODE.DistanceFunction(t, y, i, j);
            switch filtername
                case 'EnKF'
                    localization = @(t, y, H) datools.tapering.gc(t, y, r, d, H);
                case 'LETKF'
                    localization = @(t, y, Hi, k) datools.tapering.gcCTilde(t, y, Hi, r, d, k);
                case 'LETPF'
                    localization = @(t, y, Hi, k) datools.tapering.gcCTilde(t, y, Hi, r, d, k);
                case 'ETKF'
                    localization = [];
                case 'SIR'
                    localization = [];
                case 'SIS_EnKF'
                    localization = [];
                case 'ETPF'
                    localization = [];
                case 'ETPF2'
                    localization = [];
                case 'RHF'
                    localization = [];
                case 'EnGMF'
                    localization = [];
            end
        end


        filter = datools.statistical.ensemble.(filtername)('Model', model, ...
            'ObservationOperator', observationOperator, ...
            'Observation', observation, ...
            'NumEnsemble', ensN, ...
            'ModelError', modelError, ...
            'AddErrorToModel', addErrorToModel, ...
            'EnsembleGenerator', ensembleGenerator, ...
            'Inflation', inflation, ...
            'Parallel', false);


        filter.setMean(model.ODEModel.Y0);
        filter.scaleAnomalies(1/10);

        mses = zeros(steps-spinup, 1);

        rmse = nan;
        ps = '';
        doFilter = true;

        for i = 1:steps
            % forecast/propogate the truth
            truth.evolve();

            % forecast/propagate the ensembles/particles
            filter.forecast();

            % create observation based on truth (H * xt)
            xt = stateToObservation.observeWithoutError(model.TimeSpan, truth.State);
            % perturb the observation and use proper index (H * xt + noise)
            y = filter.ObservationOperator.observeWithError(model.TimeSpan, xt);

            % update the Y
            filter.Observation.updateY(y);

            % do the Analysis
            filter.analysis();

            % find the best "Estimate"
            xa = filter.BestEstimate;
            xaensembles = filter.Ensemble;

            if i > spinup
                % capture the rmse
                mses(i - spinup) = mean((xa - xt).^2);
                rmse = sqrt(mean(mses(1:(i - spinup))));
                rmse

                %rmstempvalsampleinner(i - spinup) = rmse;

            end
        end

    end


end
