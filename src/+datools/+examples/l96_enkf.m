
clear;
close all;

% time steps
Deltat =  0.05;

% Time Stepping Methods (Use ode45 or write your own)
solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 1);
solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 1);

% Define ODE 
natureODE = otp.lorenz96.presets.Canonical;
nature0 = randn(natureODE.NumVars, 1);
natureODE.TimeSpan = [0, Deltat];

modelODE = otp.lorenz96.presets.Canonical;
modelODE.TimeSpan = [0, Deltat];

% Propogate
[tt, yy] = ode45(natureODE.Rhs.F, [0 10], nature0);
natureODE.Y0 = yy(end, :).';

% initialize model
model  = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% Observation Model
naturetomodel = datools.observation.Linear(numel(nature0), 'H', speye(natureODE.NumVars));

%observeindicies = 1:natureODE.NumVars;
observeindicies = 1:1:natureODE.NumVars;

nobsvars = numel(observeindicies);

R = (1/1)*speye(nobsvars);

obserrormodel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));
%obserrormodel = datools.error.Tent;
observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indices', observeindicies);


% We make the assumption that there is no model error
modelerror = datools.error.Error;

ensembleGenerator = @(N) randn(natureODE.NumVars, N);

% ensNs = 100:100:500;
% infs = 1.01:.01:1.05;

ensNs = [25 50 75 100];
infs = [1.01 1.02 1.05 1.10];
rejs = [1.01 1.02 1.03 1.05];

% variables for which you need the rank histogram plot
histvar = 1:1:nobsvars;
rhplotval = inf*ones(numel(ensNs), numel(infs));

rmses = inf*ones(numel(ensNs), numel(infs));

maxallowerr = 10;

mm = min(rmses(:));

if  mm >= maxallowerr
    mm = 0;
end

% imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
% axis square; title('ETKF'); colormap('hot');
% xlabel('Ensemble Size'); ylabel('Inflation');

runsleft = find(rmses == inf);

f1 = figure; f2 = figure; f3 = figure; f4 = figure;

for runn = runsleft.'
    [ensNi, infi] = ind2sub([numel(ensNs), numel(infs)], runn);
    
    fprintf('N: %d, inf: %.3f\n', ensNs(ensNi), infs(infi));
    
    ns = 1;
    sE = zeros(ns, 1);
    
    inflationAll = infs(infi);
    ensN = ensNs(ensNi);
    
    for sample = 1:ns
        % Set rng for standard experiments
        rng(17 + sample - 1);
        
        inflation = inflationAll;
        
        % No localization

        r = 3;
        d = @(t, y, i, j) modelODE.DistanceFunction(t, y, i, j);
        localization = [];
        
        
        
        
        %localization = @(t, y, H) datools.tapering.gc(t, y, r, d, H);
        %$localization = @(t, y, Hi, k) datools.tapering.gcCTilde(t, y, Hi, r, d, k);
        %localization = @(t, y, Hi, k) datools.tapering.cutoffCTilde(t, y, r, d, Hi, k);
        
        enkf = datools.statistical.ensemble.ETKF(model, ...
            'Observation', observation, ...
            'NumEnsemble', ensN, ...
            'ModelError', modelerror, ...
            'EnsembleGenerator', ensembleGenerator, ...
            'Inflation', inflation, ...
            'Localization', localization, ...
            'Parallel', false, ...
            'RankHistogram', histvar, ...
            'Rejuvenation', 0.1);
        
        enkf.setMean(natureODE.Y0);
        enkf.scaleAnomalies(1/10);
        
        % define steps and spinups
        spinup = 500;
        times = 11*spinup;
        
        mses = zeros(times - spinup, 1);
        
        rmse = nan;
        
        ps = '';
        
        do_enkf = true;
        
        rmstempval = NaN * ones(1, times-spinup);
        
        for i = 1:times
            % forecast
            nature.evolve();
            
            if do_enkf
                enkf.forecast();
            end
            
            
            % observe
            xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
            y = enkf.Observation.observeWithError(model.TimeSpan(1), xt);
            
            % Rank histogram (if needed)
            datools.utils.stat.RH(enkf, xt);
            
            % analysis
            
            % try
            if do_enkf
                enkf.analysis(R, y);
            end
            %catch
            %    do_enkf = false;
            %end
            
            xa = enkf.BestEstimate;
            
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
    
    [xs, pval, rhplotval(ensNi, infi)] = datools.utils.stat.KLDiv(enkf.RankValue(1,1:end-1),...
        (1/ensN) * ones(1, ensN+1));
    
    mm = min(rmses(:));
    mm = 0;
    
    if  mm >= maxallowerr
        mm = 0;
    end
    
    
    figure(f1);
    subplot(numel(infs), numel(ensNs), runn);
    hold all;
    z = enkf.RankValue(1,1:end-1);
    z = z/sum(z);
    NN = numel(z);
    z = NN*z;
    bar(xs, z);
    plot(xs, pval, '-*r');
    set(gca,'XTick',[xs(1) xs(end)]);
    set(gca,'XTickLabel',[1, ensN+1]);
    xlabel('bins'); 
    drawnow;
    
    figure(f2);
    imagesc(ensNs.', infs.', flipud(rmses.')); caxis([0, 1]); colorbar; set(gca,'YDir','normal');
    axis square; title('EnKF'); colormap('pink');
    xlabel('Ensemble Size'); ylabel('Inflation');
    ytics = max(infs) - min(infs);
    ytics = min(infs):ytics/(length(infs) - 1):max(infs);
    set(gca,'YTick', ytics);
    set(gca,'YTickLabel', fliplr(infs));
    drawnow;
    
    figure(f3);
    imagesc(ensNs, infs, flipud(rhplotval.')); caxis([-0.09 0.09]); colorbar; set(gca, 'YDir', 'normal');
    axis square; title('KLDiv'); colormap('summer');
    xlabel('Ensemble Size'); ylabel('Inflation');
    set(gca,'YTick', ytics);
    set(gca,'YTickLabel', fliplr(infs));
    drawnow;
    
    figure(f4);
    subplot( numel(infs), numel(ensNs), runn);
    plot(spinup+1:1:times, rmstempval);
    xlim([spinup+1 times]); ylim([0 1]);
    set(gca, 'XTick', [spinup+1 times])
    set(gca, 'XTickLabel', [spinup+1 times])
    han=axes(f4,'visible','off');
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Value');
    xlabel(han,'Time Step');
    title(han,'RMSE');
    drawnow;
    
%     imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
%     axis square; title('EnKF'); colormap('pink');
%     xlabel('Ensemble Size'); ylabel('Inflation');
%     drawnow;
end

return;
