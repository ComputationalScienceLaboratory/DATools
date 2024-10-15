clear all; close all; clc;

%% User Inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this to run Lorenz96 experiments
modelname = 'Lorenz96';
% uncomment the filter you want to run
filtername = 'EnKF';
%filtername = 'ETKF';
%filtername = 'ETPF';
%filtername = 'SIR';
%filtername = 'LETKF';
%filtername = 'RHF';

%datools.statistical.ensemble.(filtername);

% observation variance
variance = 1;

% create an aray of ensemble
ensNs = [50, 75, 100, 200];

% create an array of inflation
infs = [1.02, 1.05, 1.10];

% create an array of rejuvination
rejs = 2 * logspace(-2, -1, 4);
rejs = round(rejs, 2);

% define steps and spinups
spinup = 5;
steps = 11 * spinup;

% localization radius
r = 4; % rename this
localize = true; % Change this if you need/not need localization

%% Remaining code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rank histogram for the "i" states
histvar = 1:1:1;
%decide the type of filter
switch filtername
    case 'EnKF'
        filtertype = 'Ensemble';
    case 'ETKF'
        filtertype = 'Ensemble';
    case 'LETKF'
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
Dt = 0.05;

% Time Stepping Methods (Use ode45 or write your own)
% solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 1);
% solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 1);
solvermodel = @(f, t, y) datools.utils.rk4ens(f, t, y, 1);
solvernature = @(f, t, y) datools.utils.rk4ens(f, t, y, 1);

% Define ODE for the truth
natureODE = otp.lorenz96.presets.Canonical;
nature0 = randn(natureODE.NumVars, 1);
natureODE.TimeSpan = [0, Dt];

% define ODE for the model
modelODE = otp.lorenz96.presets.Canonical;
modelODE.TimeSpan = [0, Dt];

% Propogate the truth
[tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
natureODE.Y0 = yy(end, :).';

% initialize model object
model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% Observation object that maps from nature to model
% define H. Need to have an option using linearized H (user must supply)
naturetomodel = datools.observation.Linear(numel(nature0), 'H', ...
    speye(natureODE.NumVars));

% observe these variables
% change the array if you want to observe lesser state variables
observeindicies = 1:1:natureODE.NumVars;
nobsvars = numel(observeindicies);

R = variance * speye(nobsvars);

% Observaton model object (Gaussian here)
obserrormodel = datools.uncertainty.Gaussian('Covariance', R);
observation = datools.observation.Indexed(model.NumVars, ...
    'Uncertainty', obserrormodel, ...
    'Indices', observeindicies);

% We make the assumption that there is no model error
modelerror = datools.uncertainty.NoUncertainty;

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

% f1 = figure;
% f2 = figure;
% f3 = figure;
% f4 = figure;

% total runs for all combinations of ensembles and inflation/rejuvenation
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
                filter = datools.filter.ensemble.(filtername)(model, ...
                    'InitialEnsemble', ensembleGenerator(ensN)/10, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Rejuvenation', rejuvenation);
            case 'ETKF'
                filter = datools.filter.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Rejuvenation', rejuvenation);
            case 'LETKF'
                filter = datools.filter.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Rejuvenation', rejuvenation);
            case 'SIR'
                filter = datools.filter.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Rejuvenation', rejuvenation);
            case 'ETPF'
                filter = datools.filter.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Rejuvenation', rejuvenation);
            case 'RHF'
                filter = datools.filter.ensemble.(filtername)(model, ...
                    'Observation', observation, ...
                    'NumEnsemble', ensN, ...
                    'ModelError', modelerror, ...
                    'EnsembleGenerator', ensembleGenerator, ...
                    'Inflation', inflation, ...
                    'Parallel', false, ...
                    'RankHistogram', histvar, ...
                    'Tail', 'Gaussian');

        end

        % filter.setMean(natureODE.Y0);
        % filter.scaleAnomalies(1/10);
        
        filter.MeanEstimate = natureODE.Y0;
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
            xt = naturetomodel.observeWithoutError(nature.ODEModel.TimeSpan(1), nature.State);
            % y = filter.Observation.observeWithError(model.TimeSpan(1), xt);
            y = observation.observeWithError(xt);
            observation.Uncertainty.Mean = y;


            % Rank histogram (if needed)
            datools.utils.stat.RH(filter, xt);

            % analysis
            % try
            if dofilter
                filter.analysis(observation);
            end
            %catch
            %    dofilter = false;
            %end

            xa = filter.MeanEstimate;

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

    %% Post Processing(Plotting)
    %     rw = numel(infs) - 1 - floor((runn - 1)/numel(ensNs));
    %     cl = runn - floor((runn - 1)/numel(ensNs)) * numel(ensNs);
    %
    %
    %     figure(f1);
    %     subplot(numel(infs), numel(ensNs), rw*numel(ensNs)+cl);
    %     hold all;
    %     z = filter.RankValue(1, 1:end-1);
    %     maxz = max(z);
    %     z = z / sum(z);
    %     NN = numel(z);
    %     z = NN * z;
    %     bar(xs, z);
    %     plot(xs, pval, '-*r');
    %     set(gca, 'XTick', [xs(1), xs(end)]);
    %     set(gca, 'XTickLabel', [1, ensN + 1]);
    %     set(gca, 'YTick', []);
    %     set(gca, 'YTickLabel', []);
    %     han = axes(f1, 'visible', 'off');
    %     han.Title.Visible = 'on';
    %     han.XLabel.Visible = 'on';
    %     han.YLabel.Visible = 'on';
    %     switch filtertype
    %         case 'Ensemble'
    %             ylabel(han, 'Inflation');
    %         case 'Particle'
    %             ylabel(han, 'Rejuvetion');
    %     end
    %     xlabel(han, 'Ensemble Size');
    %     title(han, 'Rank Histogram');
    %     drawnow;
    %
    %
    %     figure(f2);
    %     switch filtertype
    %         case 'Ensemble'
    %             imagesc(ensNs, infs, rmses.');
    %             caxis([0, 1]);
    %             colorbar;
    %             set(gca, 'YDir', 'normal');
    %             axis square;
    %             title('Rmse HeatMap');
    %             colormap('pink');
    %             xlabel('Ensemble Size');
    %             ylabel('Inflation')
    %             set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
    %             set(gca, 'XTickLabel', ensNs);
    %             set(gca, 'YTick', linspace(infs(1), infs(end), size(infs, 2)));
    %             set(gca, 'YTickLabel', infs);
    %             drawnow;
    %         case 'Particle'
    %             imagesc(ensNs, rejs, rmses.');
    %             caxis([0, 1]);
    %             colorbar;
    %             set(gca, 'YDir', 'normal');
    %             axis square;
    %             title('Rmse HeatMap');
    %             colormap('pink');
    %             xlabel('Ensemble Size');
    %             ylabel('Rejuvetion');
    %             set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
    %             set(gca, 'XTickLabel', ensNs);
    %             set(gca, 'YTick', linspace(rejs(1), rejs(end), size(rejs, 2)));
    %             set(gca, 'YTickLabel', rejs);
    %             drawnow;
    %     end
    %
    %     bn = bone;
    %     pk = flipud(pink);
    %     figure(f3);
    %     map1 = bn;
    %     map1 = map1(51:2:end-1, :);
    %     map2 = pk;
    %     map = [map1; map2(2:2:end-50, :)];
    %     switch filtertype
    %         case 'Ensemble'
    %             imagesc(ensNs, infs, rhplotval.');
    %             caxis([-0.1, 0.1]);
    %             colorbar;
    %             set(gca, 'YDir', 'normal');
    %             set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
    %             set(gca, 'XTickLabel', ensNs);
    %             set(gca, 'YTick', linspace(infs(1), infs(end), size(infs, 2)));
    %             set(gca, 'YTickLabel', infs);
    %             axis square;
    %             title('KLDiv');
    %             colormap(map);
    %             xlabel('Ensemble Size');
    %             ylabel('Inflation');
    %             drawnow;
    %         case 'Particle'
    %             imagesc(ensNs, rejs, rhplotval.');
    %             caxis([-0.1, 0.1]);
    %             colorbar;
    %             set(gca, 'YDir', 'normal');
    %             set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
    %             set(gca, 'XTickLabel', ensNs);
    %             set(gca, 'YTick', linspace(rejs(1), rejs(end), size(rejs, 2)));
    %             set(gca, 'YTickLabel', rejs);
    %             axis square;
    %             title('KLDiv');
    %             colormap(map);
    %             xlabel('Ensemble Size');
    %             ylabel('Rejuvetion');
    %             drawnow;
    %     end
    %
    %
    %     figure(f4);
    %     subplot(numel(infs), numel(ensNs), rw*numel(ensNs)+cl);
    %     plot(spinup+1:1:times, rmstempval);
    %     xlim([spinup + 1, times]);
    %     ylim([0, 1]);
    %     set(gca, 'XTick', [spinup + 1, times])
    %     set(gca, 'XTickLabel', [spinup + 1, times])
    %     han = axes(f4, 'visible', 'off');
    %     han.Title.Visible = 'on';
    %     han.XLabel.Visible = 'on';
    %     han.YLabel.Visible = 'on';
    %     ylabel(han, 'Value');
    %     xlabel(han, 'Time Step');
    %     title(han, 'RMSE');
    %     drawnow;

end

filename = strcat(modelname, '_', filtername, '.mat');
filepath = strcat(pwd, '\+datools\+examples\+sandu\', filename);
save(filepath);
datools.utils.plotexperiments(filepath);
return;