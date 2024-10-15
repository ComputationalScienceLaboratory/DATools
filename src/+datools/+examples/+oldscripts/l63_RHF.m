clear all; close all; clc;
figure;
drawnow;

% time steps
Delta_t = 0.12;

% Time Stepping Methods
solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 50);
solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 50);

% ODE
natureODE = otp.lorenz63.presets.Canonical;
nature0 = randn(natureODE.NumVars, 1);
natureODE.TimeSpan = [0, Delta_t];

modelODE = otp.lorenz63.presets.Canonical;
modelODE.TimeSpan = [0, Delta_t];

% Propogate
[tt, yy] = ode45(natureODE.RHS.F, [0, 10], nature0);
natureODE.Y0 = yy(end, :).';

model = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

naturetomodel = datools.observation.Linear(numel(nature0), 'H', ...
    speye(natureODE.NumVars));

%observeindicies = 1:natureODE.NumVars;
observeindicies = 1:1:natureODE.NumVars;

nobsvars = numel(observeindicies);

R = (1 / 1) * speye(nobsvars);

obserrormodel = datools.uncertainty.Gaussian('CovarianceSqrt', sqrtm(R));
%obserrormodel = datools.error.Tent;
observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indices', observeindicies);

% We make the assumption that there is no model error
modelerror = datools.error.Error;

ensembleGenerator = @(N) randn(natureODE.NumVars, N);

ensNs = 5:5:50;
infs = 1.01:.01:1.10;
serveindicies = 1:1:natureODE.NumVars;
rmses = inf * ones(numel(ensNs), numel(infs));

maxallowerr = 100;

mm = min(rmses(:));

if mm >= maxallowerr
    mm = 0;
end

imagesc(ensNs, infs, rmses.');
caxis([mm, 5]);
colorbar;
set(gca, 'YDir', 'normal');
axis square; title('RHF_l63'); colormap('hot');
xlabel('Ensemble Size');
ylabel('Inflation');

runsleft = find(rmses == inf);

for runn = runsleft.'
    [ensNi, infi] = ind2sub([numel(ensNs), numel(infs)], runn);

    fprintf('N: %d, inf: %.3f\n', ensNs(ensNi), infs(infi));

    ns = 1;
    sE = zeros(ns, 1);

    inflationAll = infs(infi);
    ensN = ensNs(ensNi);

    for sample = 1:ns
        % Set rng for standard experiments
        rng(17+sample-1);

        inflation = inflationAll;

        rhf = datools.statistical.ensemble.RHF(model, ...
            'Observation', observation, ...
            'NumEnsemble', ensN, ...
            'ModelError', modelerror, ...
            'EnsembleGenerator', ensembleGenerator, ...
            'Inflation', inflation, ...
            'Parallel', false, ...
            'Tail', 'Flat');

        %         rhf.setMean(natureODE.Y0);
        %         rhf.scaleAnomalies(1/10);

        spinup = 100;
        times = 11 * spinup;

        mses = zeros(times - spinup, 1);

        rmse = nan;

        ps = '';

        do_enkf = true;

        for i = 1:times
            % forecast
            nature.evolve();

            if do_enkf
                rhf.forecast();
            end


            % observe
            xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
            y = rhf.Observation.observeWithError(model.TimeSpan(1), xt);

            % try RH
            %datools.utils.stat.RH(enkf, xt);

            % analysis

            % try
            if do_enkf
                rhf.analysis(R, y);
            end
            %catch
            %    do_enkf = false;
            %end

            xa = rhf.BestEstimate;

            err = xt - xa;

            if i > spinup

                mses(i - spinup) = mean((xa - xt).^2);
                rmse = sqrt(mean(mses(1:(i - spinup))));

                if rmse > maxallowerr || isnan(rmse) || mses(i - spinup) > 2 * maxallowerr
                    do_enkf = false;
                end
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
    rmse
    resE = mean(sE);

    if isnan(resE)
        resE = 1000;
    end

    rmses(ensNi, infi) = resE;

    mm = min(rmses(:));
    mm = 0;

    if mm >= maxallowerr
        mm = 0;
    end

    imagesc(ensNs, infs, rmses.');
    caxis([mm, 5]);
    colorbar;
    set(gca, 'YDir', 'normal');
    axis square;
    title('RHF');
    colormap('pink');
    xlabel('Ensemble Size');
    ylabel('Inflation');
    drawnow;
end

saveas(gcf, 'plot.png')
