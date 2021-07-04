
clear;
close all;
figure;
drawnow;

Deltat =  0.11;

solvermodel = @(f, t, y) datools.utils.rk4(f, t, y, 22);
solvernature = @(f, t, y) datools.utils.rk4(f, t, y, 22);

natureODE = otp.lorenz96.presets.Canonical;
nature0 = randn(natureODE.NumVars, 1);
natureODE.TimeSpan = [0, Deltat];

modelODE = otp.lorenz96.presets.Canonical;
modelODE.TimeSpan = [0, Deltat];

[tt, yy] = ode45(natureODE.Rhs.F, [0 10], nature0);
nature0  = yy(end, :).';
natureODE.Y0 = nature0;

model  = datools.Model('Solver', solvermodel, 'ODEModel', modelODE);
nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

naturetomodel = datools.observation.Linear(numel(nature0), 'H', speye(natureODE.NumVars));

observeindicies = 2:2:natureODE.NumVars;

nobsvars = numel(observeindicies);

R = (8/1)*speye(nobsvars);

obserrormodel = datools.error.Gaussian('CovarianceSqrt', sqrtm(R));
%obserrormodel = datools.error.Tent;
observation = datools.observation.Indexed(model.NumVars, ...
    'ErrorModel', obserrormodel, ...
    'Indices', observeindicies);


% We make the assumption that there is no model error
modelerror = datools.error.Error;

ensembleGenerator = @(x) randn(natureODE.NumVars, x);

ensNs = 20:5:30;
alphass = 0.0:0.1:0.5;

rmses = inf*ones(numel(ensNs), numel(alphass));

maxallowerr = 16;

mm = min(rmses(:));

if  mm >= maxallowerr
    mm = 0;
end

runsleft = find(rmses == inf);

for runn = runsleft.'
    [ensNi, infi] = ind2sub([numel(ensNs), numel(alphass)], runn);
    
    fprintf('N: %d, inf: %.3f\n', ensNs(ensNi), alphass(infi));
    
    ns = 1;
    sE = zeros(ns, 1);
    
    alphassAll = alphass(infi);
    ensN = ensNs(ensNi);
    
    for sample = 1:ns
        % Set rng for standard experiments
        rng(17 + sample - 1);
        
        inflation = 1;
        
        % No localization
        %localization = [];
        
        %localization= @(t, y, H) datools.tapering.gc(t, y, r, d, H);
        r = 3;
        d = @(t, y, i, j) modelODE.DistanceFunction(t, y, i, j);
        localization = @(t, y, Hi, k) datools.tapering.gcCTilde(t, y, Hi, r, d, k);
        %localization = @(t, y, Hi, k) datools.tapering.cutoffCTilde(t, y, r, d, Hi, k);
        
        % required for LETPF
        localizationEnsembleDistance = @(~, ~, inds, k) spdiags([zeros(k - 1, 1); 1; zeros(numel(inds) - k + 1, 1)], 0, numel(inds), numel(inds));
        
        fnum = 1;
        
        filters = cell(fnum, 1);
        %alphas = ones(fnum*2, 1)/(2*fnum);
        alphas = ones(fnum, 1)/(fnum);
        
        for i = 1:fnum
                
            letkf = datools.statistical.ensemble.LETKF(model, ...
                'Observation', observation, ...
                'NumEnsemble', ensN, ...
                'ModelError', modelerror, ...
                'EnsembleGenerator', ensembleGenerator, ...
                'Inflation', inflation^(1/fnum), ...
                'Localization', localization, ...
                'Rejuvenation', 0, ...
                'Parallel', false);
            
            
            
            letpf = datools.statistical.ensemble.LETPF(model, ...
                'Observation', observation, ...
                'NumEnsemble', ensN, ...
                'ModelError', modelerror, ...
                'EnsembleGenerator', ensembleGenerator, ...
                'Inflation', inflation, ...
                'Localization', localization, ...
                'LocalizationEnsembleDistance', localizationEnsembleDistance, ...
                'Rejuvenation', 0, ...
                'Parallel', false);

            filters{2*i - 1} = letpf;
            filters{2*i}     = letkf;
            alphas(2*i - 1) = alphassAll;
            alphas(2*i)     = 1 - alphassAll;
                    
        end
        
        alphas = alphas/sum(alphas);
        
        
        hybrid = datools.statistical.ensemble.Hybrid(model, ...
            'Observation', observation, ...
            'NumEnsemble', ensN, ...
            'ModelError', modelerror, ...
            'EnsembleGenerator', ensembleGenerator, ...
            'Inflation', inflation, ...
            'Localization', localization, ...
            'Filters', filters, ...
            'Alphas', alphas, ...
            'Rejuvenation',  (0.2^2), ...
            'Parallel', false);
        
        
        
        natureODE.Y0 = nature0;
        nature = datools.Model('Solver', solvernature, 'ODEModel', natureODE);

        
        hybrid.setMean(natureODE.Y0);
        hybrid.scaleAnomalies(1/10);
        
        spinup = 200;
        times = 11*spinup;
        
        mses = zeros(times - spinup, 1);
        
        rmse = nan;
        
        ps = '';
        
        do_enkf = true;
        
        for i = 1:times
            % forecast
            
            nature.evolve();
            
            if do_enkf
                hybrid.forecast();
            end
            
            
            % observe
            xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
            y = hybrid.Observation.observeWithError(model.TimeSpan(1), xt);
            
            % analysis
            
            % try
            if do_enkf
                hybrid.analysis(R, y);
            end
            %catch
            %    do_enkf = false;
            %end
            
            xa = hybrid.BestEstimate;
            
            err = xt - xa;
            
            if i > spinup
                
                mses(i - spinup) = mean((xa - xt).^2);
                rmse = sqrt(mean(mses(1:(i - spinup))));
                
                if rmse > maxallowerr || isnan(rmse) || mses(i - spinup) > 200*maxallowerr
                    do_enkf = false;
                end
                
                if mod(i, 10) == 0
                    fprintf('step %d %.5f\n', i, rmse);
                end
                
            else
                
                fprintf('%d|', i)
                if mod(i, 10) == 0
                    fprintf('\n');
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
    
    imagesc(ensNs, alphass, rmses.'); caxis([1, 2.5]); colorbar; set(gca,'YDir','normal');
    axis square; title('EnKF'); colormap('pink');
    xlabel('Ensemble Size'); ylabel('\alpha');
    drawnow;
end

return;
