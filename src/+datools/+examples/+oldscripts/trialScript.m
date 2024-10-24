% trial runner script
clc;
clear;
close all;

%% Preliminaries

% ensemble, inflation and rejuvenation
ensNs = [5, 15, 25, 50];
infs = [1.01, 1.02, 1.05, 1.10];
rejs = round(logspace(-1.5, -0.25, 7), 2);

% integration time-step
deltaT = 0.12;

% filtername
filterName = 'EnKF';

%filtertype = 'Particle';
filterType = 'Ensemble';

% model
modelName = 'lorenz63';

% define steps and spinups
spinup = 500;
steps = 11 * spinup;

% Set number of samples for different initialization
numSamples = 1;

% use localization(or not) and radius
localize = false;
r = 4;

% observation error covariance
variance = 8;

modelError = datools.uncertainty.Gaussian; % we ignore model error (for now)
addErrorToModel = false; % default is false (for now)


%% decide the type of filter
switch filterName
    case 'EnKF'
        filtertype = 'Ensemble';
    case 'ETKF'
        filtertype = 'Ensemble';
    case 'LETKF'
        filtertype = 'Ensemble';
    case 'ETPF'
        filtertype = 'Particle';
    case 'ETPF2'
        filtertype = 'Particle';
    case 'LETPF'
        filtertype = 'Particle';
    case 'SIR'
        filtertype = 'Particle';
    case 'SIS_EnKF'
        filtertype = 'Particle';
    case 'RHF'
        filtertype = 'Ensemble';
    case 'EnGMF'
        filtertype = 'Ensemble';
end
fprintf('Filtername = %s, Observation Error Variance = %.2f, Runs = %d, spinups = %d\n', ...
    filterName, variance, steps, spinup);


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
odeModel.TimeSpan = [0, deltaT];

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

    switch filterType
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

    rankValueSample = zeros(numel(histVar), ensN+1, numSamples);
    rmsTempValSample = nan * ones(steps-spinup, numSamples);

    for sample = 1:numSamples
        % Set rng for standard experiments
        rng(17+sample-1);

        rmsTempValSampleInner = nan * ones(steps-spinup, 1);

        switch filterType
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
            switch filterName
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


        filter = datools.statistical.ensemble.(filterName)('Model', model, ...
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
            rankValueSample(:, :, sample) = filter.RankValue(:, 1:end-1);

            error = xa - xt;

            if i > spinup
                % capture the rmse
                mses(i - spinup) = mean((error).^2);
                rmse = sqrt(mean(mses(1:(i - spinup))));

                rmsTempValSampleInner(i - spinup) = rmse;

            end

        end
        sE(sample) = rmse;
        rmsTempValSample(:, sample) = rmsTempValSampleInner;
    end

    rankvalue = mean(rankValueSample, 3);
    rmstempval = mean(rmsTempValSample, 2);

    resE = mean(sE);

    switch filterType
        case 'Ensemble'
            rmses(ensNi, infi) = resE;

            [xs, pval, rhplotval(ensNi, infi)] = datools.utils.stat.KLDiv(rankvalue, ...
                (1 / (ensN+1))*ones(1, ensN+1));
            %rhplotval(ensNi, infi) = totalklval;
        case 'Particle'
            rmses(ensNi, reji) = resE;

            [xs, pval, rhplotval(ensNi, reji)] = datools.utils.stat.KLDiv(rankvalue, ...
                (1 / ensN)*ones(1, ensN+1));
    end

    % update all the variables for plotting
    rankvalmatrix{runn} = rankvalue;
    xvalmatrix{runn} = xs;
    polyvalmatrix{runn} = pval;
    rmsvalmatrix{runn} = rmstempval;    

end

filename = strcat(modelName, '_', filterName, '.mat');
%filepath = strcat(pwd, '\+datools\+examples\+sandu\', filename);
filepath = fullfile(pwd, '+datools', '+examples', '+v2', filename);
save(filepath);
datools.utils.plotexperiments(filepath);
