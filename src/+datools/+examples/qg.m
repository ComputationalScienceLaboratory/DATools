
clear;
close all;
figure;
drawnow;

Deltat =  0.05;

solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 1);
solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 1);

natureODE = otp.qg.presets.Canonical('Size', 'medium');
nature0 = randn(natureODE.NumVars, 1);
% natureODE.Y0 = otp.qg.presets.PopovMouIliescuSandu.relaxprolong(mic.Y0, 'medium');
natureODE.TimeSpan = [0, Deltat];

modelODE = otp.qg.presets.Canonical('Size', 'medium');
modelODE.TimeSpan = [0, Deltat];

[tt, yy] = ode45(natureODE.Rhs.F, [0 10], nature0);
natureODE.Y0 = yy(end, :).';

model  = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

naturetomodel = datools.observation.Linear(numel(nature0), 'H', speye(natureODE.NumVars));

observeindicies = 1:natureODE.NumVars;
% observeindicies = 4:4:natureODE.NumVars;

nobsvars = numel(observeindicies);

R = (1/1)*speye(nobsvars);

obserrormodel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));
%obserrormodel = datools.error.Tent;
observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indicies', observeindicies);


% We make the assumption that there is no model error
modelerror = datools.error.Error;

ensembleGenerator = @(x) randn(natureODE.NumVars, x);

ensNs = 5:5:50;
infs = 1.05:.05:1.4;

rmses = inf*ones(numel(ensNs), numel(infs));

maxallowerr = 2;

mm = min(rmses(:));

if  mm >= maxallowerr
    mm = 0;
end

imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
axis square; title('ETKF'); colormap('hot');
xlabel('Ensemble Size'); ylabel('Inflation');

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
        rng(17 + sample - 1);
        
        inflation = inflationAll;
        
        % No localization
        r = 4;
        d = @(t, y, i, j) modelODE.DistanceFunction(t, y, i, j);
        localization = [];
        
        localization= @(t, y, H) datools.tapering.gc(t, y, r, d, H);
        
        %localization = @(t, y, Hi, k) datools.tapering.gcCTilde(t, y, r, d, Hi, k);
        %localization = @(t, y, Hi, k) datools.tapering.cutoffCTilde(t, y, r, d, Hi, k);
        
        enkf = datools.statistical.ensemble.EnKF(model, ...
            'Observation', observation, ...
            'NumEnsemble', ensN, ...
            'ModelError', modelerror, ...
            'EnsembleGenerator', ensembleGenerator, ...
            'Inflation', inflation, ...
            'Localization', localization, ...
            'Parallel', false);
        
        enkf.setMean(nature0);
        enkf.scaleAnomalies(1/10);
        
        spinup = 100;
        times = 11*spinup;
        
        mses = zeros(times - spinup, 1);
        
        rmse = nan;
        
        ps = '';
        
        do_enkf = true;
        
        for i = 1:times
            % forecast
            
            nature.evolve();
            
            if do_enkf
                enkf.forecast();
            end
            
            
            % observe
            xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
            y = enkf.Observation.observeWithError(model.TimeSpan(1), xt);
            
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
                
                if rmse > maxallowerr || isnan(rmse) || mses(i - spinup) > 2*maxallowerr
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
    
    resE = mean(sE);
    
    if isnan(resE)
        resE = 1000;
    end
    
    rmses(ensNi, infi) = resE;
    
    mm = min(rmses(:));
    
    if  mm >= maxallowerr
        mm = 0;
    end
    
    imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
    axis square; title('EnKF'); colormap('hot');
    xlabel('Ensemble Size'); ylabel('Inflation');
    drawnow;
end

return;