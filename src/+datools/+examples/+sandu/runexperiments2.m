function runexperiments2(user)
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

% numberof samples for averaging runs
numsamples = user.ns;

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
    case 'ETPF2'
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
        nature0 = randn(natureODE.NumVars, 1); % natureODE.NumVars are the number of variables for the model
        natureODE.TimeSpan = [0, Dt];
        
        modelODE = otp.lorenz63.presets.Canonical;
        modelODE.TimeSpan = [0, Dt];
        
        % Propogate the model
        [tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
        natureODE.Y0 = yy(end, :).';
        
        xt = natureODE.Y0;
        
        % This can be used to generate ensemble if needed
        ensembleGenerator = @(N) xt + 0.1*randn(natureODE.NumVars, N);
        
        % observe these variables
        % change the array if you want to observe lesser state variables
        observeindicies = 1:1:natureODE.NumVars;
        nobsvars = numel(observeindicies);   
    case 'Lorenz96'
        natureODE = otp.lorenz96.presets.Canonical;
        nature0 = randn(natureODE.NumVars, 1); % natureODE.NumVars are the number of variables for the model
        natureODE.TimeSpan = [0, Dt];
        
        modelODE = otp.lorenz96.presets.Canonical;
        modelODE.TimeSpan = [0, Dt];
        
        % Propogate the model
        [tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
        natureODE.Y0 = yy(end, :).';
        
        xt = natureODE.Y0;
        
        % This can be used to generate ensemble if needed
        ensembleGenerator = @(N) xt + 0.1*randn(natureODE.NumVars, N);
        
                % observe these variables
        % change the array if you want to observe lesser state variables
        observeindicies = 1:1:natureODE.NumVars;
        nobsvars = numel(observeindicies);
    case 'QG'
        natureODE = otp.qg.presets.Canonical('Size', [63, 127]);
        nature0 = natureODE.Y0;
        natureODE.TimeSpan = [0, Dt];
        
        modelODE = otp.qg.presets.Canonical('Size', [63, 127]);
        modelODE.TimeSpan = [0, Dt];
        
        xt = natureODE.Y0;
        
        load('qgtrajectory.mat')
        y = y.';
        Nt = size(y, 2);
        
        ensembleGenerator = @(N) xt + 0.1*y(:, randperm(Nt, N));
        
        % observe these variables
        % change the array if you want to observe lesser state variables
        nobsvars = 150;
        observeindicies = round(linspace(1, natureODE.NumVars, nobsvars));
end

% nature0 = randn(natureODE.NumVars, 1); % natureODE.NumVars are the number of variables for the model
% natureODE.TimeSpan = [0, Dt];


%modelODE.TimeSpan = [0, Dt];

% Propogate the model
% [tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
% natureODE.Y0 = yy(end, :).';

% initialize model
model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% Observation Model
naturetomodel = datools.observation.Linear(numel(nature0), 'H', ...
    speye(natureODE.NumVars));



R = variance * speye(nobsvars);

% Observaton model (Gaussian here)
obserrormodel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));

observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indices', observeindicies);

% We make the assumption that there is no model error
modelerror = datools.error.Error;

serveindicies = 1:1:natureODE.NumVars;
switch(filtertype)
    case 'Ensemble'
        rmses = inf * ones(numel(ensNs), numel(infs));
        rhplotval = inf * ones(numel(ensNs), numel(infs));
    case 'Particle'
        rmses = inf * ones(numel(ensNs), numel(rejs));
        rhplotval = inf * ones(numel(ensNs), numel(rejs));
end
% rmses = inf * ones(numel(ensNs), numel(infs));
% rhplotval = inf * ones(numel(ensNs), numel(infs));
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

    ns = numsamples;
    sE = zeros(ns, 1);
    rankvaluesample = zeros(histvar, ensN+1, numsamples);
    rmstempvalsample = nan * ones(steps-spinup, numsamples);

    parfor sample = 1:ns
        %fprintf('N: %d | inf: %.3f | sample = %d\n', ensNs(ensNi), infs(infi), sample);
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
                inflation = 0;
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
                case 'ETPF'
                    localization = [];
                case 'ETPF2'
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
            case 'ETPF2'
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

        %rmstempval = nan * ones(1, steps-spinup);

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
            rankvaluesample(:,:,sample) = filter.RankValue(1, 1:end-1); 

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

                %rmstempvalsample(i - spinup,sample) = rmse;
                
                rmstempvalsampleinner(i - spinup) = rmse;

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
        
        rmstempvalsample(:,sample) = rmstempvalsampleinner;
        
    end
    
    rankvalue = mean(rankvaluesample,3);
    rmstempval = mean(rmstempvalsample,2);

    resE = mean(sE);

    if isnan(resE)
        resE = 1000;
    end

    switch filtertype
        case 'Ensemble'
            rmses(ensNi, infi) = resE;

            [xs, pval, rhplotval(ensNi, infi)] = datools.utils.stat.KLDiv(rankvalue, ...
                (1 / ensN)*ones(1, ensN+1));
        case 'Particle'
            rmses(ensNi, reji) = resE;

            [xs, pval, rhplotval(ensNi, reji)] = datools.utils.stat.KLDiv(rankvalue, ...
                (1 / ensN)*ones(1, ensN+1));
    end

    mm = min(rmses(:));
    mm = 0;

    if mm >= maxallowerr
        mm = 0;
    end

    % update all the variables for plotting
    rankvalmatrix{runn} = rankvalue;
    xvalmatrix{runn} = xs;
    polyvalmatrix{runn} = pval;
    rmsvalmatrix{runn} = rmstempval;

    totalruns = totalruns + 1;
end

filename = strcat(modelname, '_', filtername, '.mat');
%filepath = strcat(pwd, '\+datools\+examples\+sandu\', filename);
filepath = fullfile(pwd, '+datools', '+examples', '+sandu', filename);
save(filepath);
datools.utils.plotexperiments(filepath);
end