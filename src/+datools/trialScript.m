% trial runner script
clc;
clear all;
close all;

%% Preliminaries

% ensemble, inflation and rejuvenation
ensNs = [5, 15, 25, 50];
infs = [1.01, 1.02, 1.05, 1.10];
rejs = round(logspace(-1.5, -0.25, 7), 2);

% integration time-step
Delta_t = 0.05;

% filtername
filtername = 'LETPF2';

filtertype = 'Particle';
%filtertype = 'Ensemble';


% model
modelName = 'lorenz96';

% define steps and spinups
spinup = 50;
steps = 11 * spinup;

% Set number of samples for different initialization
numSamples = 1;

% use localization(or not) and radius
localize = false;
r = 4;

% observation error covariance
variance = 8;

modelError = datools.error.Gaussian; % we ignore model error (for now)
addErrorToModel = false; % default is false (for now)


%% decide the type of filter
% switch filtername
%     case 'EnKF'
%         filtertype = 'Ensemble';
%     case 'ETKF'
%         filtertype = 'Ensemble';
%     case 'LETKF'
%         filtertype = 'Ensemble';
%     case 'ETPF'
%         filtertype = 'Particle';
%     case 'ETPF2'
%         filtertype = 'Particle';
%     case 'LETPF'
%         filtertype = 'Particle';
%     case 'SIR'
%         filtertype = 'Particle';
%     case 'SIS_EnKF'
%         filtertype = 'Particle';
%     case 'RHF'
%         filtertype = 'Ensemble';
%     case 'EnGMF'
%         filtertype = 'Ensemble';
% end
fprintf('Filtername = %s, Observation Error Variance = %.2f, Runs = %d, spinups = %d\n', ...
    filtername, variance, steps, spinup);


%% Preliminaries variables for plotting
histVar = 1:1:1;
measureRankHist = 'Truth';

%plotting parameters
rankhistogramplotindex = 1:2:numel(ensNs);
rmseplotindex = 1:2:numel(ensNs);
rmseheatmapplotindex = 1:2:numel(ensNs);
kldivergenceplotindex = 1:2:numel(ensNs);

% save the necessary variables
rmsvalmatrix = {};
rankvalmatrix = {};
xvalmatrix = {};
polyvalmatrix = {};


%% define the ODE Model
odeOTP = otp.(modelName).presets.Canonical;
odeModel = datools.ODEModel('OTPObject', odeOTP);
odeModel.TimeSpan = [0, Delta_t];

% define ODE solvers (time integrators) for truth and model
solverModel = @(f, t, y) ode45(f, t, y);
solverTruth = @(f, t, y) ode45(f, t, y);

%define the truth and model ODE
model = datools.Model('Solver', solverModel, 'ODEModel', odeModel, 'SynthError', modelError, ...
    'AddError', addErrorToModel);
truth = datools.Model('Solver', solverTruth, 'ODEModel', odeModel, 'SynthError', modelError, ...
    'AddError', addErrorToModel);


%% define observation and model error

% observe these variables
observeIndicies = 1:2:3; % change accordingly

numStateObserved = numel(observeIndicies);

R = variance * speye(numStateObserved);

% Observaton error model (Gaussian here)
obsErrorModel = datools.error.Gaussian('Covariance', R);

% Observation Operator & data (use this for assimilation)
observation = datools.observation.operator.IndexedObservation(model.NumVars, ...
    'ErrorModel', obsErrorModel, ...
    'NumObs', numStateObserved, ...
    'Indices', observeIndicies);

% Observation
% observation = datools.observation.StateObservation('Covariance', R, 'Noise', obsErrorModel, ...
%     'NumObs', numStateObserved);

% Observation model to map from state space to observation space (use this for truth observation)
stateToObservation = datools.observation.operator.LinearObservation(truth.NumVars, ...
    'H', speye(truth.NumVars));

% create random ensemble generator
ensembleGenerator = @(N) randn(model.NumVars, N);


%% iterate 
rmses = inf * ones(numel(ensNs), numel(infs));

runsLeft = find(rmses == inf);



for runn = runsLeft.'

    switch filtertype
        case 'Ensemble'
            [ensNi, infi] = ind2sub([numel(ensNs), numel(infs)], runn);
            reji = 0;
            inflationAll = infs(infi);
            rejAll = 0;
            ensN = ensNs(ensNi);
        case 'Particle'
            [ensNi, reji] = ind2sub([numel(ensNs), numel(rejs)], runn);
            infi = 0;
            rejAll = rejs(reji);
            inflationAll = 0;
            ensN = ensNs(ensNi);
    end

    %[ensNi, infi] = ind2sub([numel(ensNs), numel(infs)], runn);

    sE = zeros(numSamples, 1);

    rankvaluesample = zeros(numel(histVar), ensN+1, numSamples);
    rmstempvalsample = nan * ones(steps-spinup, numSamples);

    for sample = 1:numSamples
        % Set rng for standard experiments
        rng(17+sample-1);

        rmstempvalsampleinner = nan * ones(steps-spinup, 1);

        switch filtertype
            case 'Ensemble'
                inflation = inflationAll;
                rejuvenation = 0;
                fprintf('N: %d | inf: %.3f | sample = %d\n', ensNs(ensNi), infs(infi), sample);
            case 'Particle'
                rejuvenation = rejAll;
                inflation = 1.05;
                fprintf('N: %d | rej: %.3f | sample = %d\n', ensNs(ensNi), rejs(reji), sample);
        end

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
            'ModelError', modelError, ...
            'AddErrorToModel', addErrorToModel, ...
            'NumEnsemble', ensN, ...
            'EnsembleGenerator', ensembleGenerator, ...
            'Inflation', inflation, ...
            'Rejuvenation', rejuvenation, ...
            'Parallel', false, ...
            'RankHistogram', histVar, ...
            'Tail', 'Gaussian');


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

            % create observation based on truth (H * xt) full space!
            xt = stateToObservation.observeWithoutError(model.TimeSpan(1), truth.State);

            % perturb the observation and use proper index (H * xt + noise)
            y = observation.observeWithError(model.TimeSpan(1), xt);

            % update the Y of the observation object
            observation.updateY(y);

            % do the analysis/assimilation
            filter.analysis(observation);

            % find the best "Estimate"
            xa = filter.BestEstimate;
            xaEnsembles = filter.Ensemble;

            % observable
            hxa = stateToObservation.observeWithoutError(model.TimeSpan(1), xaEnsembles); % H(x_analysis)
            Hxa = observation.observeWithError(model.TimeSpan(1), hxa); % H(x_analysis) + noise

            % Rank histogram (if needed)
            datools.utils.stat.RH(filter, observation, xt, y, Hxa, measureRankHist);
            rankvaluesample(:, :, sample) = filter.RankValue(:, 1:end-1);

            error = xa - xt;

            if i > spinup
                % capture the rmse
                mses(i - spinup) = mean((error).^2);
                rmse = sqrt(mean(mses(1:(i - spinup))));
                rmse

                %rmstempvalsampleinner(i - spinup) = rmse;

            end

        end

    end


end
