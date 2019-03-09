
addpath('/home/apopov/Dropbox/Projects/ode-test-problems/src');

csl.datools.presetmodels.mlqgsoEn

inflation = sqrt(1.6);
%inflation = sqrt((ensN - 1)/(ensN - 3));

%inflation = 1;

%radius = 15.0;
radius = 16;
localization = @(t, y, H) csl.datools.tapering.gc(t, y, radius, distfn, H);
%localization = @(t, y, H) csl.datools.tapering.gauss(t, y, radius, distfn, H);


% enkfstart = csl.datools.statistical.ensemble.DEnKF(model{end}, ...
%     'Observation', observation, ...
%     'NumEnsemble', ensN*numel(model), ...
%     'ModelError', modelerror, ...
%     'EnsembleGenerator', ensembleGenerator, ...
%     'Inflation', inflation, ...
%     'Localization', localization, ...
%     'Parallel', true);

% starttimes = 10;
% 
% for i = 1:starttimes
%     
%     % forecast
%     nature.evolve();
%     enkfstart.forecast();
% 
%     % observe
%     xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
%     y = enkfstart.Observation.observeWithError(model{1}.TimeSpan(1), xt);
%     
%     % analysis
%     enkfstart.analysis(R, y);
% 
% end
% 
% 
% ensembleGenerator = @(x) enkfstart.Ensemble(:, randperm(ensN*numel(model), x));


%model = {model{2}}

enkf = csl.datools.statistical.ensemble.MLDEnKF(model, ...
    'Observation', observation, ...
    'NumEnsemble', ensN, ...
    'ModelError', modelerror, ...
    'EnsembleGenerator', ensembleGenerator, ...
    'Inflation', inflation, ...
    'Localization', localization, ...
    'Parallel', true);

spinup = 500;
times = 11*spinup;

spinup = 5;
times = 1000;



mses = zeros(times - spinup, 1);

rmse = nan;

ps = '';

for i = 1:times
    
    % forecast
    nature.evolve();
    enkf.forecast();

    % observe
    xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
    y = enkf.Observation.observeWithError(model{1}.TimeSpan(1), xt);
    
    % analysis
    enkf.analysis(R, y);
    
    xa = enkf.BestEstimate;
    
    subplot(2, 1, 1);
    imagesc(reshape(xt, 127, 127));
    axis square; colorbar;
    subplot(2, 1, 2);
    imagesc(reshape(xa, 127, 127));
    axis square; colorbar;
    drawnow;
    
    if i > spinup
        mses(i - spinup) = mean((xa - xt).^2);
        rmse = sqrt(mean(mses(1:(i - spinup))));
    end
    
    for kk = 1:numel(ps)
        fprintf('\b');
    end
    
    ps = sprintf('step: %d, rmse: %.5f\n', i, rmse);
    
    fprintf(ps);
end
