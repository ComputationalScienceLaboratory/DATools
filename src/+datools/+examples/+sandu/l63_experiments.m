clear all; close all;

%% Preprocessing(User Inputs)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this to run Lorenz 63 experiments

% uncomment the filter you want to run
filtername = 'EnKF';
%filtername = 'ETKF';
%filtername = 'ETPF';
%filtername = 'SIR';
%filtername = 'RHF';

% oservation variance
variance = 1;

% create an aray of ensemble
ensNs = [25, 50, 75, 100];

% create an array of inflation
infs = [1.01, 1.02, 1.05, 1.10];

% create an array of rejuvination
rejs = 2 * logspace(-2, -1, 4);
rejs = round(rejs, 2);

% define steps and spinups
spinup = 50;
times = 11 * spinup;

% Mention the location and name to save (change this)(uncomment the last line)
savdir = '/home/abhinab93/Documents/experiments/Lorenz63/ETPF/l63ETPF.mat';

%% Remaining code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    filtername, variance, times, spinup);
% time steps
Delta_t = 0.12;

% Time Stepping Methods (Use ode45 or write your own)
solvermodel = @(f, t, y) ode45(f, t, y);
solvernature = @(f, t, y) ode45(f, t, y);

% Define ODE (using OTP)
natureODE = otp.lorenz63.presets.Canonical;
nature0 = randn(natureODE.NumVars, 1); % natureODE.NumVars are the number of variables for the model
natureODE.TimeSpan = [0, Delta_t];

modelODE = otp.lorenz63.presets.Canonical;
modelODE.TimeSpan = [0, Delta_t];

% Propogate the model
[tt, yy] = ode45(natureODE.Rhs.F, [0, 10], nature0);
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

R = variance * (1 / 1) * speye(nobsvars);

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

maxallowerr = 10;

mm = min(rmses(:));

if mm >= maxallowerr
    mm = 0;
end

runsleft = find(rmses == inf);

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

% total runs for all combinations of ensembles and inflation/rejuvenation
for runn = runsleft.'
    switch filtertype
        case 'Ensemble'
            [ensNi, infi] = ind2sub([numel(ensNs), numel(infs)], runn);
            fprintf('N: %d, inf: %.3f\n', ensNs(ensNi), infs(infi));
            inflationAll = infs(infi);
        case 'Particle'
            [ensNi, reji] = ind2sub([numel(ensNs), numel(rejs)], runn);
            fprintf('N: %d, rej: %.3f\n', ensNs(ensNi), rejs(reji));
            rejAll = rejs(reji);
    end

    ns = 1;
    sE = zeros(ns, 1);

    ensN = ensNs(ensNi);

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

        end


        filter.setMean(natureODE.Y0);
        filter.scaleAnomalies(1/10);

        mses = zeros(times - spinup, 1);

        rmse = nan;
        ps = '';
        do_filter = true;

        rmstempval = NaN * ones(1, times-spinup);

        % Assimilation steps
        for i = 1:times
            % forecast/evolve the model
            nature.evolve();

            if do_filter
                filter.forecast();
            end

            % observe
            xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
            y = filter.Observation.observeWithError(model.TimeSpan(1), xt);

            % Rank histogram (if needed)
            datools.utils.stat.RH(filter, xt);

            % analysis

            % try
            if do_filter
                filter.analysis(R, y);
            end
            %catch
            %    do_enkf = false;
            %end

            xa = filter.BestEstimate;

            err = xt - xa;

            if i > spinup
                % capture the rmse
                mses(i - spinup) = mean((xa - xt).^2);
                rmse = sqrt(mean(mses(1:(i - spinup))));

                rmstempval(i - spinup) = rmse;

                %                 if rmse > maxallowerr || isnan(rmse) || mses(i - spinup) > 2*maxallowerr
                %                     do_enkf = false;
                %                 end
            end

            if ~do_filter
                break;
            end

        end
        hold off;

        if isnan(rmse)
            rmse = 1000;
        end

        if ~do_filter
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

    %% Post Processing(Plotting)
    rw = numel(infs) - 1 - floor((runn - 1)/numel(ensNs));
    cl = runn - floor((runn - 1)/numel(ensNs)) * numel(ensNs);


    figure(f1);
    subplot(numel(infs), numel(ensNs), rw*numel(ensNs)+cl);
    hold all;
    z = filter.RankValue(1, 1:end-1);
    maxz = max(z);
    z = z / sum(z);
    NN = numel(z);
    z = NN * z;
    bar(xs, z);
    plot(xs, pval, '-*r');
    set(gca, 'XTick', [xs(1), xs(end)]);
    set(gca, 'XTickLabel', [1, ensN + 1]);
    set(gca, 'YTick', []);
    set(gca, 'YTickLabel', []);
    han = axes(f1, 'visible', 'off');
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    switch filtertype
        case 'Ensemble'
            ylabel(han, 'Inflation');
        case 'Particle'
            ylabel(han, 'Rejuvetion');
    end
    xlabel(han, 'Ensemble Size');
    title(han, 'Rank Histogram');
    drawnow;


    figure(f2);
    switch filtertype
        case 'Ensemble'
            imagesc(ensNs, infs, rmses.');
            caxis([0, 1]);
            colorbar;
            set(gca, 'YDir', 'normal');
            axis square;
            title('Rmse HeatMap');
            colormap('pink');
            xlabel('Ensemble Size');
            ylabel('Inflation')
            set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
            set(gca, 'XTickLabel', ensNs);
            set(gca, 'YTick', linspace(infs(1), infs(end), size(infs, 2)));
            set(gca, 'YTickLabel', infs);
            drawnow;
        case 'Particle'
            imagesc(ensNs, rejs, rmses.');
            caxis([0, 1]);
            colorbar;
            set(gca, 'YDir', 'normal');
            axis square;
            title('Rmse HeatMap');
            colormap('pink');
            xlabel('Ensemble Size');
            ylabel('Rejuvetion');
            set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
            set(gca, 'XTickLabel', ensNs);
            set(gca, 'YTick', linspace(rejs(1), rejs(end), size(rejs, 2)));
            set(gca, 'YTickLabel', rejs);
            drawnow;
    end

    bn = bone;
    pk = flipud(pink);
    figure(f3);
    map1 = bn;
    map1 = map1(51:2:end-1, :);
    map2 = pk;
    map = [map1; map2(2:2:end-50, :)];
    switch filtertype
        case 'Ensemble'
            imagesc(ensNs, infs, rhplotval.');
            caxis([-0.1, 0.1]);
            colorbar;
            set(gca, 'YDir', 'normal');
            set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
            set(gca, 'XTickLabel', ensNs);
            set(gca, 'YTick', linspace(infs(1), infs(end), size(infs, 2)));
            set(gca, 'YTickLabel', infs);
            axis square;
            title('KLDiv');
            colormap(map);
            xlabel('Ensemble Size');
            ylabel('Inflation');
            drawnow;
        case 'Particle'
            imagesc(ensNs, rejs, rhplotval.');
            caxis([-0.1, 0.1]);
            colorbar;
            set(gca, 'YDir', 'normal');
            set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
            set(gca, 'XTickLabel', ensNs);
            set(gca, 'YTick', linspace(rejs(1), rejs(end), size(rejs, 2)));
            set(gca, 'YTickLabel', rejs);
            axis square;
            title('KLDiv');
            colormap(map);
            xlabel('Ensemble Size');
            ylabel('Rejuvetion');
            drawnow;
    end


    figure(f4);
    subplot(numel(infs), numel(ensNs), rw*numel(ensNs)+cl);
    plot(spinup+1:1:times, rmstempval);
    xlim([spinup + 1, times]);
    ylim([0, 1]);
    set(gca, 'XTick', [spinup + 1, times])
    set(gca, 'XTickLabel', [spinup + 1, times])
    han = axes(f4, 'visible', 'off');
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabel(han, 'Value');
    xlabel(han, 'Time Step');
    title(han, 'RMSE');
    drawnow;

end


%save(fullfile(savdir));
return;