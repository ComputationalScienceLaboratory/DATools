clear; close all;
% figure;
% drawnow;

% time steps
Delta_t = 0.12;

% Time Stepping Methods (Use ode45 or write your own)
% solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 50);
% solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 50);
solvermodel = @(f, t, y) ode45(f, t, y);
solvernature = @(f, t, y) ode45(f, t, y);

% Define ODE 
natureODE = otp.lorenz63.presets.Canonical;
nature0 = randn(natureODE.NumVars, 1);
natureODE.TimeSpan = [0, Delta_t];

modelODE = otp.lorenz63.presets.Canonical;
modelODE.TimeSpan = [0, Delta_t];

% Propogate
[tt, yy] = ode45(natureODE.Rhs.F, [0 10], nature0);
natureODE.Y0 = yy(end, :).';

% initialize model
model  = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

% Observation Model
naturetomodel = datools.observation.Linear(numel(nature0), 'H',...
    speye(natureODE.NumVars));


% observe these variables
observeindicies = 1:1:natureODE.NumVars;

nobsvars = numel(observeindicies);

R = (1/1)*speye(nobsvars);

% Observaton model (Gaussian here)
obserrormodel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));
%obserrormodel = datools.error.Tent;
observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indices', observeindicies);

% We make the assumption that there is no model error
modelerror = datools.error.Error;

ensembleGenerator = @(N) randn(natureODE.NumVars, N);

% number of ensembles and inflation
% ensNs = 10:5:25;
% infs = 1.01:0.01:1.04;

%ensNs = [5 15 25 50];
ensNs = [50 100 150 200];
infs = [1.01 1.02 1.05 1.10];
rejs = [0.10 0.12 0.15 0.20];

% variables for which you need the rank histogram plot
histvar = 1:1:3;
serveindicies = 1:1:natureODE.NumVars;
rmses = inf*ones(numel(ensNs), numel(infs));
rhplotval = inf*ones(numel(ensNs), numel(infs));

maxallowerr = 10;

mm = min(rmses(:));

if  mm >= maxallowerr
    mm = 0;
end

% imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
% axis square; title('EnKF_l63'); colormap('hot');
% xlabel('Ensemble Size'); ylabel('Inflation');

runsleft = find(rmses == inf);

f1 = figure; f2 = figure; f3 = figure; f4 = figure;

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
        rng(17 + sample - 1);
        
        inflation = inflationAll;
        rejuvenation = rejAll;
        
        % define the statistical/variational model here
        enkf = datools.statistical.ensemble.SIR(model, ...
            'Observation', observation, ...
            'NumEnsemble', ensN, ...
            'ModelError', modelerror, ...
            'EnsembleGenerator', ensembleGenerator, ...
            'Inflation', inflation, ...
            'Parallel', false, ...
            'RankHistogram', histvar, ...
            'Rejuvenation', rejuvenation);
        
        enkf.setMean(natureODE.Y0);
        enkf.scaleAnomalies(1/10);
        
        % define steps and spinups
        spinup = 100;
        times = 11*spinup;
        
        mses = zeros(times - spinup, 1);
        
        %         imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
        %         axis square; title('EnKF'); colormap('pink');
        %         xlabel('Ensemble Size'); ylabel('Inflation');
        
        rmse = nan;
        ps = '';
        do_enkf = true;
        
        rmstempval = NaN * ones(1, times-spinup);
        
        % Assimilation
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
        hold off;
        
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
    maxz = max(z);
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
    %imagesc(ensNs.', infs.', flipud(rmses.')); caxis([0, 1]); colorbar; set(gca,'YDir','normal');
    imagesc(ensNs.', rejs.', flipud(rmses.')); caxis([0, 1]); colorbar; set(gca,'YDir','normal');
    axis square; title('ETPF'); colormap('pink');
    xlabel('Ensemble Size'); ylabel('Inflation');
    %ytics = max(infs) - min(infs);
    ytics = max(rejs) - min(rejs);
    %ytics = min(infs):ytics/(length(infs) - 1):max(infs);
    ytics = min(rejs):ytics/(length(rejs) - 1):max(rejs);
    set(gca,'YTick', ytics);
    %set(gca,'YTickLabel', fliplr(infs));
    set(gca,'YTickLabel', fliplr(rejs));
    drawnow;
    
    figure(f3);
    %imagesc(ensNs.', infs.', flipud(rhplotval.')); caxis([-0.09 0.09]); colorbar; set(gca, 'YDir', 'normal');
    imagesc(ensNs, rejs, flipud(rhplotval.')); caxis([-0.09 0.09]); colorbar; set(gca, 'YDir', 'normal');
    axis square; title('KLDiv'); colormap('summer');
    xlabel('Ensemble Size'); ylabel('Inflation');
    set(gca,'YTick', ytics);
    %set(gca,'YTickLabel', fliplr(infs));
    set(gca,'YTickLabel', fliplr(rejs));
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
    
end
% step = 0;
% Rank histogram
% figure;
% for i = 1: length(enkf.RankValue(:,1))
%     subplot(1,3,i);
%     bar(enkf.RankValue(i,1:end-1));
% end

% Kldiv + poly approx (for the variable being observed)
% figure;
% for i = 1:length(histvar)
%     
% end


return;
