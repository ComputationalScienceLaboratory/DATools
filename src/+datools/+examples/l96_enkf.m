
clear;
close all;
figure;
drawnow;

Deltat =  0.05;

solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 1);
solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 1);

natureODE = otp.lorenz96.presets.Canonical;
nature0 = randn(natureODE.NumVars, 1);
natureODE.TimeSpan = [0, Deltat];

modelODE = otp.lorenz96.presets.Canonical;
modelODE.TimeSpan = [0, Deltat];

[tt, yy] = ode45(natureODE.Rhs.F, [0 10], nature0);
natureODE.Y0 = yy(end, :).';

model  = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

naturetomodel = datools.observation.Linear(numel(nature0), 'H', speye(natureODE.NumVars));

%observeindicies = 1:natureODE.NumVars;
observeindicies = 1:1:natureODE.NumVars;

nobsvars = numel(observeindicies);

R = (1/1)*speye(nobsvars);

obserrormodel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));
%obserrormodel = datools.error.Tent;
observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indicies', observeindicies);


% We make the assumption that there is no model error
modelerror = datools.error.Error;

ensembleGenerator = @(N) randn(natureODE.NumVars, N);

ensNs = 100:100:500;
infs = 1.01:.01:1.05;
infs = 0.1:0.1:0.7;

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

        %r = 3;
        %d = @(t, y, i, j) modelODE.DistanceFunction(t, y, i, j);
        localization = [];
        
        %localization = @(t, y, H) datools.tapering.gc(t, y, r, d, H);
        

        %localization = @(t, y, Hi, k) datools.tapering.gcCTilde(t, y, Hi, r, d, k);
        %localization = @(t, y, Hi, k) datools.tapering.cutoffCTilde(t, y, r, d, Hi, k);
        
        enkf = datools.statistical.ensemble.SIR(model, ...
            'Observation', observation, ...
            'NumEnsemble', ensN, ...
            'ModelError', modelerror, ...
            'EnsembleGenerator', ensembleGenerator, ...
            'Inflation', inflation, ...
            'Localization', localization, ...
            'Parallel', false);
        
        enkf.setMean(natureODE.Y0);
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
    mm = 0;
    
    if  mm >= maxallowerr
        mm = 0;
    end
    
    imagesc(ensNs, infs, rmses.'); caxis([mm, 1]); colorbar; set(gca,'YDir','normal');
    axis square; title('EnKF'); colormap('pink');
    xlabel('Ensemble Size'); ylabel('Inflation');
    drawnow;
end

return;
