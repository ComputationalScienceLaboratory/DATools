clear;
close all;
clc;

% time steps
dt = 0.05;
filtername = 'EnKF';

% Time Stepping Methods (Use ode45 or write your own)
solvermodel = @(f, t, y) datools.utils.rk4ens(f, t, y, 1);
solvernature = @(f, t, y) datools.utils.rk4ens(f, t, y, 1);

% Define ODE for truth 
otpNature = otp.lorenz63.presets.Canonical;
% natureODE = otp.lorenz96.presets.Canonical;
natureODE = datools.ODEModel('OTPObject', otpNature);
nature0 = randn(natureODE.NumVars, 1);
natureODE.TimeSpan = [0, dt];

% Define ODE for the model
otpModel = otp.lorenz63.presets.Canonical;
% modelODE = otp.lorenz96.presets.Canonical;
modelODE = datools.ODEModel('OTPObject', otpModel);
modelODE.TimeSpan = [0, dt];

% Propogate
% [tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
[tt, yy] = ode45(natureODE.F, [0, 10], nature0);
natureODE.Y0 = yy(end, :).';

% initialize model
model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% Observation Model
naturetomodel = @(x) x;

%observeindicies = 1:natureODE.NumVars;
observeindicies = 1:1:natureODE.NumVars;

nobsvars = numel(observeindicies);

R = (1 / 1) * speye(nobsvars);
% R = (1 / 1) * eye(nobsvars); % for EnGMF


obserrormodel = datools.uncertainty.Gaussian('Covariance', R);
observation = datools.observation.Indexed(model.NumVars, ...
    'Uncertainty', obserrormodel, ...
    'Indices', observeindicies);

natureobserrormodel = datools.uncertainty.Gaussian('Covariance', R);
natureobs = datools.observation.Indexed(model.NumVars, ...
    'Uncertainty', natureobserrormodel, ...
    'Indices', observeindicies);


% We make the assumption that there is no model error
modelerror = datools.uncertainty.NoUncertainty;

ensembleGenerator = @(N) randn(natureODE.NumVars, N);

% ensNs = 100:100:500;
% infs = 1.01:.01:1.05;

ensNs = [15, 25, 50, 100];
infs = [1.01, 1.02, 1.05, 1.10];
rejs = 2 * logspace(-2, -1, 4);
rejs = round(rejs, 2);

% variables for which you need the rank histogram plot
histvar = 1:1:nobsvars;
rhplotval = inf * ones(numel(ensNs), numel(infs));

rmses = inf * ones(numel(ensNs), numel(infs));

maxallowerr = 10;

mm = min(rmses(:));

if mm >= maxallowerr
    mm = 0;
end

% imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
% axis square; title('ETKF'); colormap('hot');
% xlabel('Ensemble Size'); ylabel('Inflation');

runsleft = find(rmses == inf);

f1 = figure;
f2 = figure;
f3 = figure;
f4 = figure;

for runn = runsleft.'
    [ensNi, infi] = ind2sub([numel(ensNs), numel(infs)], runn);
    [ensNi, reji] = ind2sub([numel(ensNs), numel(rejs)], runn);

    fprintf('N: %d, inf: %.3f\n', ensNs(ensNi), infs(infi));

    ns = 1;
    sE = zeros(ns, 1);

    inflationAll = infs(infi);
    rejAll = rejs(reji);
    ensN = ensNs(ensNi);

    for sample = 1:ns
        % Set rng for standard experiments
        rng(17+sample-1);

        inflation = inflationAll;
        rejuvenation = rejAll;

        % No localization

        r = 4;
        d = @(y, i, j) modelODE.DistanceFunction(0, y, i, j);
        

        localization = [];
        % localization = @(y, H) datools.tapering.bloc.gc(y, r, d, H);
        % localization = @(y, H, k) datools.tapering.rloc.gc(y, r, d, H, k);
        % localization = @(t, y, Hi, k) datools.tapering.gcCTilde(t, y, Hi, r, d, k);
        % localization = @(t, y, Hi, k) datools.tapering.cutoffCTilde(t, y, r, d, Hi, k);


        filter = datools.filter.ensemble.(filtername)(model, ...
            'InitialEnsemble', ensembleGenerator(ensN)/10, ...
            'Inflation', inflation, ...
            'Localization', localization, ...
            'Parallel', false, ...
            'RankHistogram', histvar, ...
            'Rejuvenation', rejuvenation);

        filter.MeanEstimate = natureODE.Y0;

        % define steps and spinups
        spinup = 500;
        times = 11 * spinup;

        mses = zeros(times - spinup, 1);

        rmse = nan;

        ps = '';

        do_enkf = true;

        rmstempval = NaN * ones(1, times-spinup);

        for i = 1:times
            % forecast
            nature.evolve();

            if do_enkf
                filter.forecast();
            end


            % observe
            xt = naturetomodel(nature.State);
            y = natureobs.observeWithError(xt);
            observation.Uncertainty.Mean = y;

            % Rank histogram (if needed)
            datools.utils.stat.RH(filter, xt, y);

            % analysis

            % try
            if do_enkf
                filter.analysis(observation);
            end
            %catch
            %    do_enkf = false;
            %end

            xa = filter.MeanEstimate;

            err = xt - xa;

            if i > spinup

                mses(i - spinup) = mean((xa - xt).^2);
                rmse = sqrt(mean(mses(1:(i - spinup))));

                rmstempval(i - spinup) = rmse;

                %                 if rmse > maxallowerr || isnan(rmse) || mses(i - spinup) > 2*maxallowerr
                %                     do_enkf = false;
                %                 end
            end


            if ~do_enkf
                break;
            end

        end

        if isnan(rmse)
            rmse = 1000;
        end

        if ~do_enkf
            sE(sample) = 1000;
        else
            sE(sample) = rmse;
        end

    end

    resE = mean(sE);

    if isnan(resE)
        resE = 1000;
    end

    rmses(ensNi, infi) = resE;

    [xs, pval, rhplotval(ensNi, infi)] = datools.utils.stat.KLDiv(filter.RankValue(1, 1:end-1), ...
        (1 / ensN)*ones(1, ensN+1));

    mm = min(rmses(:));
    mm = 0;

    if mm >= maxallowerr
        mm = 0;
    end

    rw = numel(infs) - 1 - floor((runn - 1)/numel(ensNs));
    cl = runn - floor((runn - 1)/numel(ensNs)) * numel(ensNs);
    figure(f1);
    set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
    set(gca, 'XTickLabel', ensNs);
    subplot(numel(infs), numel(ensNs), rw*numel(ensNs)+cl);
    hold all;
    z = filter.RankValue(1, 1:end-1);
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
    %han.XTick.Visible = 'on';
    %han.XTickLabel.Visible = 'on';
    %set(han, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs,2)));
    %set(han, 'XTickLabel', ensNs);
    ylabel(han, 'Inflation');
    xlabel(han, 'Ensemble Size');
    title(han, 'Rank Histogram');
    drawnow;

    figure(f2);
    imagesc(ensNs, infs, rmses.');
    caxis([0, 1]);
    colorbar;
    set(gca, 'YDir', 'normal');
    set(gca, 'XTick', linspace(ensNs(1), ensNs(end), size(ensNs, 2)));
    set(gca, 'XTickLabel', ensNs);
    set(gca, 'YTick', linspace(infs(1), infs(end), size(infs, 2)));
    set(gca, 'YTickLabel', infs);
    axis square;
    title('Rmse HeatMap');
    colormap('pink');
    xlabel('Ensemble Size');
    ylabel('Inflation');
    drawnow;

    figure(f3);
    map = bone;
    map = map(1:2:end-1, :);
    pt = flipud(pink);
    map = [map; pt(2:2:end, :)];
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
    xlabel(han, 'Time Steps');
    title(han, 'RMSE');
    drawnow;

    %     imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
    %     axis square; title('EnKF'); colormap('pink');
    %     xlabel('Ensemble Size'); ylabel('Inflation');
    %     drawnow;
end
% savdir = '/home/abhinab93/Documents/experiments/Lorenz96/EnKFr4/l96enkf.mat';
% save(fullfile(savdir));
return;
