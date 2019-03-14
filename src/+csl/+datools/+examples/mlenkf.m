clear;
% Set rng for standard experiments
rng(20);

addpath('../../ode-test-problems/src');

csl.datools.presetmodels.mlqgsoEn

%inflation = sqrt(1.08);
inflation = sqrt((ensN - 1)/(ensN - 3));

%inflation = 1;

%radius = 25.0;
radius = 5;
%localization = @(t, y, H) csl.datools.tapering.gc(t, y, radius, distfn, H);
localization{2} = @(t, y, H) csl.datools.tapering.gauss(t, y, radius, distfn, H);

radiuslow = 5;
localization{1} = @(t, y, H) csl.datools.tapering.gauss(t, y, radiuslow, distfn, H);




%model = {model{end}};

enkf = csl.datools.statistical.ensemble.MLEnKF(model, ...
    'Observation', observation, ...
    'NumEnsemble', ensN, ...
    'ModelError', modelerror, ...
    'EnsembleGenerator', ensembleGenerator, ...
    'Inflation', inflation, ...
    'Localization', localization, ...
    'Parallel', true);



enkfC = csl.datools.statistical.ensemble.POEnKF(modelC, ...
    'Observation', observation, ...
    'NumEnsemble', ensN, ...
    'ModelError', modelerror, ...
    'EnsembleGenerator', ensembleGenerator, ...
    'Inflation', inflation, ...
    'Localization', localization{2}, ...
    'Parallel', true);

enkfC.Ensemble = enkf.Ensembles{end};

%spinup = 500;
%times = 11*spinup;

spinup = 25;
times = 250;



msesML = zeros(times - spinup, 1);
msesC = zeros(times - spinup, 1);

rmseML = nan;
rmseC = nan;

ps = '';

for i = 1:times
    
    % forecast
    nature.evolve();
    enkf.forecast();
    enkfC.forecast();

    % observe
    xt = naturetomodel.observeWithoutError(nature.TimeSpan(1), nature.State);
    y = enkf.Observation.observeWithError(model{1}.TimeSpan(1), xt);
    
    % analysis
    enkf.analysis(R, y);
    enkfC.analysis(R, y);
    
    xa = enkf.BestEstimate;
    xaC = enkfC.BestEstimate;
    
    subplot(2, 1, 1);
    imagesc(reshape(xt, 127, 127));
    axis square; colorbar;
    subplot(2, 2, 3);
    imagesc(reshape(xa, 127, 127));
    axis square; colorbar;
    title('ML');
    subplot(2, 2, 4);
    imagesc(reshape(xaC, 127, 127));
    axis square; colorbar;
    title('Control');
    drawnow;
    
    if i > spinup
        msesML(i - spinup) = mean((xa - xt).^2);
        rmseML = sqrt(mean(msesML(1:(i - spinup))));
        
        msesC(i - spinup) = mean((xaC - xt).^2);
        rmseC = sqrt(mean(msesC(1:(i - spinup))));
    end
    
    for kk = 1:numel(ps)
        fprintf('\b');
    end
    
    ps = sprintf('step: %d, rmseML: %.5f, rmseC: %.5f\n', i, rmseML, rmseC);
    
    fprintf(ps);
end
