function runexperiments(user)
% model name
modelname = user.modelname;

% filtername
filtername = user.filtername;

% oservation variance
variance = user.variance;

% Number of Ensembles in the model run
ensNs = user.ensNs;

% Inflation (if required)
infs = user.infs;

% Rejuvenation(if required)
rejs = user.rejs;

% Steps and Spinups
spinup = user.spinups;
steps = user.steps;

% localization
localize = user.localize;

%% Remaining code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rank histogram for the 1st state only (for now)
histvar = 1:1:1;

% decide the type of filter
switch filtername
    case 'EnKF'
        filtertype = 'Ensemble';
    case 'ETKF'
        filtertype = 'Ensemble';
    case 'ETPF'
        filtertype = 'Particle';
    case 'SIR'
        filtertype = 'Particle';
    case 'RHF'
        filtertype = 'Ensemble';
end
fprintf('Filtername = %s, Observation Variance = %.2f, Runs = %d, spinups = %d\n', ...
    filtername, variance, steps, spinup);

% time steps
Dt = user.Dt;

% Time Stepping Methods (Use ode45 or write your own)
odesolver = user.odesolver;
switch odesolver
    case 'ode45'
        solvermodel = @(f, t, y) ode45(f, t, y);
        solvernature = @(f, t, y) ode45(f, t, y);
    case 'RK4'
        solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 1);
        solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 1);
end

% Model
switch modelname
    case 'Lorenz63'
        natureODE = otp.lorenz63.presets.Canonical;
        modelODE = otp.lorenz63.presets.Canonical;
    case 'Lorenz96'
        natureODE = otp.lorenz96.presets.Canonical;
        modelODE = otp.lorenz96.presets.Canonical;
end

nature0 = randn(natureODE.NumVars, 1); % natureODE.NumVars are the number of variables for the model
natureODE.TimeSpan = [0, Dt];


modelODE.TimeSpan = [0, Dt];

% Propogate the model
[tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
natureODE.Y0 = yy(end, :).';

% initialize model
model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% Observation Model
naturetomodel = datools.observation.Linear(numel(nature0), 'H', ...
    speye(natureODE.NumVars));

% observe these variables
% change the array if you want to observe lesser state variables
observeindicies = 1:1:natureODE.NumVars;
nobsvars = numel(observeindicies);

R = variance * speye(nobsvars);

% Observaton model (Gaussian here)
obserrormodel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));

observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indices', observeindicies);

% We make the assumption that there is no model error
modelerror = datools.error.Error;

% This can be used to generate ensemble if needed
ensembleGenerator = @(N) randn(natureODE.NumVars, N);

serveindicies = 1:1:natureODE.NumVars;
rmses = inf * ones(numel(ensNs), numel(infs));
rhplotval = inf * ones(numel(ensNs), numel(infs));
totalruns = 0;

% save the necessary variables
rmsvalmatrix = {};
rankvalmatrix = {};
xvalmatrix = {};
polyvalmatrix = {};

maxallowerr = 10;

mm = min(rmses(:));

if mm >= maxallowerr
    mm = 0;
end

runsleft = find(rmses == inf);


for runn = runsleft.'
    switch filtertype
        case 'Ensemble'
            [ensNi, infi] = ind2sub([numel(ensNs), numel(infs)], runn);
            fprintf('N: %d, inf: %.3f\n', ensNs(ensNi), infs(infi));
            inflationAll = infs(infi);
            ensN = ensNs(ensNi);
        case 'Particle'
            [ensNi, reji] = ind2sub([numel(ensNs), numel(rejs)], runn);
            fprintf('N: %d, rej: %.3f\n', ensNs(ensNi), rejs(reji));
            rejAll = rejs(reji);
            ensN = ensNs(ensNi);
    end

    ns = 1;
    sE = zeros(ns, 1);

    for sample = 1:ns
        % Set rng for standard experiments
        rng(17+sample-1);

        switch filtertype
            case 'Ensemble'
                inflation = inflationAll;
                rejuvenation = 0;
            case 'Particle'
                rejuvenation = rejAll;
                inflation = 0;
        end

        if (localize)
            d = @(t, y, i, j) modelODE.DistanceFunction(t, y, i, j);
            switch filtername
                case 'EnKF'
                    localization = @(t, y, H) datools.tapering.gc(t, y, r, d, H);
                case 'LETKF'
                    localization = @(t, y, Hi, k) datools.tapering.gcCTilde(t, y, Hi, r, d, k);
                case 'ETKF'
                    localization = [];
                case 'SIR'
                    localization = [];
                case 'ETPF'
                    localization = [];
                case 'RHF'
                    localization = [];
            end
        end

        switch filtername
            case 'EnKF'
                filter = datools.statistical.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar);
            case 'ETKF'
                filter = datools.statistical.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar);
            case 'ETPF'
                filter = datools.statistical.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Rejuvenation', rejuvenation);
            case 'SIR'
                filter = datools.statistical.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Rejuvenation', rejuvenation);
            case 'RHF'
                filter = datools.statistical.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Tail', 'Gaussian');
            case 'LETKF'
                filter = datools.statistical.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Rejuvenation', rejuvenation);

        end

        filter.setMean(natureODE.Y0);
        filter.scaleAnomalies(1/10);

        mses = zeros(steps-spinup, 1);

        rmse = nan;
        ps = '';
        dofilter = true;

        rmstempval = NaN * ones(1, steps-spinup);

        % Assimilation steps
        for i = 1:steps
            % forecast/evolve the model
            nature.evolve();

            if dofilter
                filter.forecast();
            end

            % observe
            xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
            y = filter.Observation.observeWithError(model.TimeSpan(1), xt);

            % Rank histogram (if needed)
            datools.utils.stat.RH(filter, xt);

            % analysis
            % try
            if dofilter
                filter.analysis(R, y);
            end
            %catch
            %    dofilter = false;
            %end

            xa = filter.BestEstimate;

            err = xt - xa;

            if i > spinup
                % capture the rmse
                mses(i - spinup) = mean((xa - xt).^2);
                rmse = sqrt(mean(mses(1:(i - spinup))));

                rmstempval(i - spinup) = rmse;

                %                 if rmse > maxallowerr || isnan(rmse) || mses(i - spinup) > 2*maxallowerr
                %                     dofilter = false;
                %                 end
            end

            if ~dofilter
                break;
            end

        end

        if isnan(rmse)
            rmse = 1000;
        end

        if ~dofilter
            sE(sample) = 1000;
        else
            sE(sample) = rmse;
        end

    end

    resE = mean(sE);

    if isnan(resE)
        resE = 1000;
    end

    switch filtertype
        case 'Ensemble'
            rmses(ensNi, infi) = resE;

            [xs, pval, rhplotval(ensNi, infi)] = datools.utils.stat.KLDiv(filter.RankValue(1, 1:end-1), ...
                (1 / ensN)*ones(1, ensN+1));
        case 'Particle'
            rmses(ensNi, reji) = resE;

            [xs, pval, rhplotval(ensNi, reji)] = datools.utils.stat.KLDiv(filter.RankValue(1, 1:end-1), ...
                (1 / ensN)*ones(1, ensN+1));
    end

    mm = min(rmses(:));
    mm = 0;

    if mm >= maxallowerr
        mm = 0;
    end

    % update all the variables for plotting
    rankvalmatrix{runn} = filter.RankValue(1, 1:end-1);
    xvalmatrix{runn} = xs;
    polyvalmatrix{runn} = pval;
    rmsvalmatrix{runn} = rmstempval;

    totalruns = totalruns + 1;
end

filename = strcat(modelname, '_', filtername, '.mat');
filepath = strcat(pwd, '\+datools\+examples\+sandu\', filename);
save(filepath);
datools.utils.plotexperiments(filepath);
end